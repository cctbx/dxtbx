from __future__ import absolute_import, division, print_function

import h5py
import libtbx
from scitbx.array_family import flex
from dxtbx.format.FormatNexus import FormatNexus
from dxtbx.model import MultiAxisGoniometer

try:
    from dxtbx.util.masking import GoniometerMaskerFactory
except ImportError:
    GoniometerMaskerFactory = False


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

        if not GoniometerMaskerFactory:
            return False

        # Get the file handle
        handle = h5py.File(image_file, "r")
        if "short_name" not in handle["/entry/instrument"].attrs:
            return False
        if handle["/entry/instrument"].attrs["short_name"].lower() not in (
            b"i03",
            b"i04",
        ):
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
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        super(FormatNexusEigerDLS16M, self).__init__(image_file, **kwargs)
        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)

    def get_mask(self, index, goniometer=None):
        mask = super(FormatNexusEigerDLS16M, self).get_mask()
        if mask is None:
            # XXX mask really shouldn't be None
            # https://jira.diamond.ac.uk/browse/SCI-8308
            mask = tuple(
                flex.bool(flex.grid(reversed(panel.get_image_size())), True)
                for panel in self.get_detector()
            )
        if self._dynamic_shadowing and self.get_scan():
            gonio_masker = self.get_goniometer_shadow_masker(goniometer=goniometer)
            scan = self.get_scan()
            detector = self.get_detector()
            shadow_mask = gonio_masker.get_mask(
                detector, scan.get_angle_from_image_index(index)
            )
            assert len(mask) == len(shadow_mask)
            for m, sm in zip(mask, shadow_mask):
                if sm is not None:
                    m &= sm
        return mask

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
