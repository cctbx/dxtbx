from __future__ import annotations

import math

import h5py

from libtbx import Auto

import dxtbx.nexus.nxmx
from dxtbx.format.FormatNXmx import FormatNXmx
from dxtbx.format.FormatNXmxDLS import FormatNXmxDLS
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
