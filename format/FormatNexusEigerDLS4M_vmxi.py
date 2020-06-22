from __future__ import absolute_import, division, print_function

import ast
import sys

import h5py

import libtbx

from dxtbx.format.FormatNexus import FormatNexus
from dxtbx.masking import GoniometerMaskerFactory
from dxtbx.model import MultiAxisGoniometer


def get_count_limit_from_meta(meta_file_name):
    with h5py.File(meta_file_name, "r") as f:

        config = f["/config"][()]
        config_data = ast.literal_eval(config.decode("utf-8"))

    return config_data["countrate_correction_count_cutoff"]


class FormatNexusEigerDLS4M_vmxi(FormatNexus):
    @staticmethod
    def understand(image_file):
        with h5py.File(image_file, "r") as handle:
            name = FormatNexusEigerDLS4M_vmxi.get_instrument_name(handle)
            if name is None or name.lower() != b"vmxi":
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

        super(FormatNexusEigerDLS4M_vmxi, self).__init__(image_file, **kwargs)
        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)

    def get_detector(self, index=None):
        # workaround for https://jira.diamond.ac.uk/browse/I03-365
        # read the count limit from the meta file - if anything goes
        # wrong, do nothing

        detector = self._detector()

        try:
            if self._image_file.endswith("_master.h5"):
                meta = self._image_file.replace("_master.h5", "_meta.h5")
            elif self._image_file.endswith(".nxs"):
                meta = self._image_file.replace(".nxs", "_meta.h5")
            limit = get_count_limit_from_meta(meta)

            assert limit > 0

            for panel in detector:
                trusted = panel.get_trusted_range()
                panel.set_trusted_range((trusted[0], limit))
        except Exception:
            pass

        return detector

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
    for arg in sys.argv[1:]:
        print(FormatNexusEigerDLS4M_vmxi.understand(arg))
