from __future__ import annotations

import math

import h5py
import nxmx

from libtbx import Auto

from dxtbx.format.FormatNXmx import FormatNXmx
from dxtbx.format.FormatNXmxDLS import FormatNXmxDLS, get_bit_depth_from_meta
from dxtbx.masking import GoniometerMaskerFactory


class FormatNXmxDLSI19_2(FormatNXmxDLS):
    """
    Read NXmx-flavour data from beamline I19-2 at Diamond Light Source.

    Include the option of dynamic shadowing of the standard I19-2 diamond anvil
    pressure cell with a 76° conical aperture.
    """

    @staticmethod
    def understand(image_file):
        """This format class applies if the instrument name contains 'I19-2'."""
        with h5py.File(image_file, swmr=True) as handle:
            name = nxmx.h5str(FormatNXmx.get_instrument_name(handle))
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
        try:
            self._bit_depth_readout = get_bit_depth_from_meta(self._meta)
        except Exception:
            self._bit_depth_readout = 16

    def get_goniometer_shadow_masker(self, goniometer=None):
        """Apply the dynamic mask for a diamond anvil cell with a 76° aperture."""
        if goniometer is None:
            goniometer = self.get_goniometer()

        return GoniometerMaskerFactory.diamond_anvil_cell(
            goniometer, cone_opening_angle=math.radians(76)
        )

    # workaround for negative rotation angles used in serial experiments
    def _start(self):
        self._static_mask = None

        self._cached_file_handle = h5py.File(self._image_file, swmr=True)
        nxmx_obj = self._get_nxmx(self._cached_file_handle)
        nxentry = nxmx_obj.entries[0]
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
            nxdata = nxmx_obj.entries[0].data[0]
            if nxdata.signal:
                data = nxdata[nxdata.signal]
            else:
                data = list(nxdata.values())[0]
            self._num_images, *_ = data.shape
