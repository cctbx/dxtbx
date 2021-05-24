"""
An implementation of the CBF image reader for Pilatus images, from the Pilatus
6M SN 100 currently on Diamond I04.
"""


import binascii
import sys

import libtbx
from scitbx.array_family import flex

from dxtbx.ext import uncompress
from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.masking import GoniometerMaskerFactory


class FormatCBFMiniPilatusDLS6MSN100(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 100 @ DLS."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

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
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)
        super().__init__(image_file, **kwargs)

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header, though for miniCBF
        there are limited options for this."""

        if "Phi" in self._cif_header_dictionary:
            phi_value = float(self._cif_header_dictionary["Phi"].split()[0])
        else:
            phi_value = 0.0

        if "Kappa" in self._cif_header_dictionary:
            kappa_value = float(self._cif_header_dictionary["Kappa"].split()[0])
        else:
            kappa_value = 0.0

        if "Omega" in self._cif_header_dictionary:
            omega_value = float(self._cif_header_dictionary["Omega"].split()[0])
        else:
            omega_value = 0.0

        phi = (1.0, 0.0, 0.0)
        kappa = (0.914, 0.279, -0.297)
        omega = (1.0, 0.0, 0.0)
        axes = flex.vec3_double((phi, kappa, omega))
        angles = flex.double((phi_value, kappa_value, omega_value))
        names = flex.std_string(("GON_PHI", "GON_KAPPA", "GON_OMEGA"))
        return self._goniometer_factory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    def get_goniometer_shadow_masker(self, goniometer=None):
        return GoniometerMaskerFactory.mini_kappa(goniometer)

    def _read_cbf_image(self):
        start_tag = binascii.unhexlify("0c1a04d5")

        with self.open_file(self._image_file, "rb") as fh:
            data = fh.read()
        data_offset = data.find(start_tag) + 4
        cbf_header = self._parse_cbf_header(
            data[: data_offset - 4].decode("ascii", "ignore")
        )

        pixel_values = uncompress(
            packed=data[data_offset : data_offset + cbf_header["size"]],
            fast=cbf_header["fast"],
            slow=cbf_header["slow"],
        )

        return pixel_values


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDLS6MSN100.understand(arg))
