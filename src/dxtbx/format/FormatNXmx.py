from __future__ import annotations

import h5py

import dxtbx.nexus
from dxtbx.format.FormatNexus import FormatNexus


class FormatNXmx(FormatNexus):
    """Read NXmx-flavour NeXus-format HDF5 data."""

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            return bool(
                [
                    entry
                    for entry in dxtbx.nexus.nxmx.find_class(handle, "NXentry")
                    if "definition" in entry
                    and dxtbx.nexus.nxmx.h5str(entry["definition"][()]) == "NXmx"
                ]
            )

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        super().__init__(image_file, **kwargs)

    def _start(self):
        self._static_mask = None

        self._cached_file_handle = h5py.File(self._image_file, swmr=True)
        nxmx = self._get_nxmx(self._cached_file_handle)
        nxentry = nxmx.entries[0]
        nxsample = nxentry.samples[0]
        nxinstrument = nxentry.instruments[0]
        nxdetector = nxinstrument.detectors[0]
        nxbeam = nxinstrument.beams[0]
        self._goniometer_model = dxtbx.nexus.get_dxtbx_goniometer(nxsample)
        self._beam_factory = dxtbx.nexus.CachedWavelengthBeamFactory(nxbeam)
        wavelength = self._beam_factory.make_beam(index=0).get_wavelength()
        self._detector_model = dxtbx.nexus.get_dxtbx_detector(nxdetector, wavelength)
        self._scan_model = dxtbx.nexus.get_dxtbx_scan(nxsample, nxdetector)
        self._static_mask = dxtbx.nexus.get_static_mask(nxdetector)
        self._bit_depth_readout = nxdetector.bit_depth_readout

        if self._scan_model:
            self._num_images = len(self._scan_model)
        else:
            nxdata = nxmx.entries[0].data[0]
            if nxdata.signal:
                data = nxdata[nxdata.signal]
            else:
                data = list(nxdata.values())[0]
            self._num_images, *_ = data.shape

    def _get_nxmx(self, fh: h5py.File):
        return dxtbx.nexus.nxmx.NXmx(fh)

    def _beam(self, index: int | None = None) -> dxtbx.model.Beam:
        return self._beam_factory.make_beam(index=index or 0)

    def get_spectrum(self, index=None):
        return self._beam_factory.make_spectrum(index=index or 0)

    def get_num_images(self) -> int:
        return self._num_images

    def get_static_mask(self, index=None, goniometer=None):
        return self._static_mask

    def get_raw_data(self, index):
        nxmx = self._get_nxmx(self._cached_file_handle)
        nxdata = nxmx.entries[0].data[0]
        nxdetector = nxmx.entries[0].instruments[0].detectors[0]
        raw_data = dxtbx.nexus.get_raw_data(nxdata, nxdetector, index)
        if self._bit_depth_readout:
            # if 32 bit then it is a signed int, I think if 8, 16 then it is
            # unsigned with the highest two values assigned as masking values
            if self._bit_depth_readout == 32:
                top = 2**31
            else:
                top = 2**self._bit_depth_readout
            for data in raw_data:
                d1d = data.as_1d()
                d1d.set_selected(d1d == top - 1, -1)
                d1d.set_selected(d1d == top - 2, -2)
        return raw_data
