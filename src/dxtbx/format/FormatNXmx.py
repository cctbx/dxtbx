from __future__ import annotations

import math

import h5py

from libtbx import Auto

import dxtbx.nexus
from dxtbx.format.FormatNexus import FormatNexus
from dxtbx.masking import GoniometerMaskerFactory


class FormatNXmx(FormatNexus):
    """Read NXmx-flavour NeXus-format HDF5 data."""

    _cached_file_handle = None

    @staticmethod
    def understand(image_file):
        """This format class currently only applies to beamline I19-2 at DLS."""
        with h5py.File(image_file, swmr=True) as handle:
            name = dxtbx.nexus.nxmx.h5str(FormatNXmx.get_instrument_name(handle))
        if name and "I19-2" in name:
            return True
        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        super().__init__(image_file, **kwargs)

    def _start(self):
        self._static_mask = None

        with h5py.File(self._image_file, swmr=True) as fh:
            nxmx = dxtbx.nexus.nxmx.NXmx(fh)
            nxsample = nxmx.entries[0].samples[0]
            nxinstrument = nxmx.entries[0].instruments[0]
            nxdetector = nxinstrument.detectors[0]
            nxbeam = nxinstrument.beams[0]

            self._goniometer_model = dxtbx.nexus.get_dxtbx_goniometer(nxsample)
            self._beam_model = dxtbx.nexus.get_dxtbx_beam(nxbeam)
            self._detector_model = dxtbx.nexus.get_dxtbx_detector(nxdetector, nxbeam)
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

    def _beam(self, index=None):
        return self._beam_model

    def get_num_images(self) -> int:
        return self._num_images

    def get_static_mask(self, index=None, goniometer=None):
        return self._static_mask

    def get_raw_data(self, index):
        if self._cached_file_handle is None:
            self._cached_file_handle = h5py.File(self._image_file, swmr=True)

        nxmx = dxtbx.nexus.nxmx.NXmx(self._cached_file_handle)
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


class FormatNXmxI19_2(FormatNXmx):
    """
    Read NXmx-flavour data from beamline I19-2 at Diamond Light Source.

    Include the option of dynamic shadowing of the standard I19-2 diamond anvil
    pressure cell with a 76° conical aperture.
    """

    @staticmethod
    def understand(image_file):
        """This format class applies if the instrument name contains 'I19-2'."""
        with h5py.File(image_file, swmr=True) as handle:
            name = dxtbx.nexus.nxmx.h5str(FormatNXmx.get_instrument_name(handle))
        if name and "I19-2" in name:
            return True
        return False

    @staticmethod
    def has_dynamic_shadowing(**kwargs):
        """Check if dynamic shadowing should be applied for a diamond anvil cell."""
        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if dynamic_shadowing in (Auto, "Auto"):
            return False
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)
        super().__init__(image_file, **kwargs)

    def get_goniometer_shadow_masker(self, goniometer=None):
        """Apply the dynamic mask for a diamond anvil cell with a 76° aperture."""
        if goniometer is None:
            goniometer = self.get_goniometer()

        return GoniometerMaskerFactory.diamond_anvil_cell(
            goniometer, cone_opening_angle=math.radians(76)
        )
