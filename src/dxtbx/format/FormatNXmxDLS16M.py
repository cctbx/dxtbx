import h5py

import libtbx

from dxtbx.format.FormatNXmxDLS import FormatNXmxDLS
from dxtbx.masking import GoniometerMaskerFactory
from dxtbx.model import MultiAxisGoniometer

# These are the instrument names that should be used according to
# https://manual.nexusformat.org/classes/applications/NXmx.html and
# https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_diffrn_source.type.html
VALID_NAMES = {
    # "long" names
    "DIAMOND BEAMLINE I03",
    "DIAMOND BEAMLINE I04",
    # "short" names
    "DLS I03",
    "DLS I04",
    # "legacy" names used until 2021/03/12
    "I03",
    "I04",
}


class FormatNXmxDLS16M(FormatNXmxDLS):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            name = FormatNXmxDLS16M.get_instrument_name(handle)
            return name and name.upper() in VALID_NAMES

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

        super().__init__(image_file, **kwargs)
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
                f"Don't understand this goniometer: {list(goniometer.get_names())}"
            )
