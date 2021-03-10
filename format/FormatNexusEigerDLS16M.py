from __future__ import absolute_import, division, print_function

import h5py

import libtbx

from dxtbx.format.FormatNexusEigerDLS import FormatNexusEigerDLS
from dxtbx.masking import GoniometerMaskerFactory
from dxtbx.model import MultiAxisGoniometer

VALID_NAMES = {
    # "long" names
    b"DIAMOND BEAMLINE I03",
    b"DIAMOND BEAMLINE I04",
    # "short" names
    b"DLS I03",
    b"DLS I04",
}

LEGACY_NAMES = {b"I03", b"I04"}


class FormatNexusEigerDLS16M(FormatNexusEigerDLS):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            name = FormatNexusEigerDLS.get_instrument_name(handle)
            if name is None:
                return False
            if name.upper() in LEGACY_NAMES:
                return True
            if name in VALID_NAMES:
                return True

        return False

    def has_dynamic_shadowing(self, **kwargs):
        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if not isinstance(self.get_goniometer(), MultiAxisGoniometer):
            # Single-axis goniometer, no goniometer shadows
            return False
        if dynamic_shadowing in (libtbx.Auto, "Auto"):
            return True
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        super(FormatNexusEigerDLS16M, self).__init__(image_file, **kwargs)
        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)

    def get_goniometer_shadow_masker(self, goniometer=None):
        if not self._dynamic_shadowing:
            return None

        if goniometer is None:
            goniometer = self.get_goniometer()

        assert goniometer is not None

        if goniometer.get_names()[1] == "chi":
            return GoniometerMaskerFactory.smargon(goniometer)
        elif goniometer.get_names()[1] == "kappa":
            return GoniometerMaskerFactory.mini_kappa(goniometer)

        else:
            raise RuntimeError(
                "Don't understand this goniometer: %s" % list(goniometer.get_names())
            )
