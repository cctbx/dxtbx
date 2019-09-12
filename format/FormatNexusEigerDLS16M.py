from __future__ import absolute_import, division, print_function

import h5py
import libtbx
from dxtbx.format.FormatNexus import FormatNexus
from dxtbx.model import MultiAxisGoniometer
from dxtbx.masking import GoniometerMaskerFactory


class FormatNexusEigerDLS16M(FormatNexus):
    @staticmethod
    def understand(image_file):
        """Check to see if this format class can understand the image file.

        Args:
          image_file (str): The file path of the image file to check.

        Returns:
          bool: Returns ``True`` if the image_file is understood by this format class,
          else returns ``False``.

        """

        # this depends on DIALS for the goniometer shadow model; if missing
        # simply return False

        # Get the file handle
        handle = h5py.File(image_file, "r")
        name = FormatNexusEigerDLS16M.get_instrument_name(handle)
        if name is None or name.lower() not in (b"i03", b"i04"):
            return False

        return True

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


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatNexusEigerDLS16M.understand(arg))
