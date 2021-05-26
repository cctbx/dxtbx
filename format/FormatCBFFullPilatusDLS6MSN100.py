import sys

import libtbx

from dxtbx.format.FormatCBFFullPilatus import FormatCBFFullPilatus
from dxtbx.masking import GoniometerMaskerFactory


class FormatCBFFullPilatusDLS6MSN100(FormatCBFFullPilatus):
    """An image reading class for full CBF format images from Pilatus
    detectors."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        # this depends on DIALS for the goniometer shadow model; if missing
        # simply return False

        header = FormatCBFFullPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 60-0100 Diamond" in header
            ):
                return True

        return False

    @staticmethod
    def has_dynamic_shadowing(**kwargs):
        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if dynamic_shadowing in (libtbx.Auto, "Auto"):
            return True
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)
        super().__init__(image_file, **kwargs)

    def get_goniometer_shadow_masker(self, goniometer=None):
        if goniometer is None:
            goniometer = self.get_goniometer()

        assert goniometer is not None
        if goniometer.get_names()[1] == "GON_CHI":
            return GoniometerMaskerFactory.smargon(goniometer)

        elif goniometer.get_names()[1] == "GON_KAPPA":
            return GoniometerMaskerFactory.mini_kappa(goniometer)

        else:
            raise ValueError(
                "Don't understand this goniometer: %s" % list(goniometer.get_names())
            )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFFullPilatus.understand(arg))
