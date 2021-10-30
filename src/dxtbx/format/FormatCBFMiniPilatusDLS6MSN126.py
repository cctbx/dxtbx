"""
An implementation of the CBF image reader for Pilatus images, from the Pilatus
6M SN 100 currently on Diamond I03.
"""


import binascii
import sys

from scitbx import matrix

from dxtbx.ext import uncompress
from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus


class FormatCBFMiniPilatusDLS6MSN126(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 126 @ DLS."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 60-0126" in header
            ):
                return True

        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""
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

        axis = matrix.col((1, 0, 0))
        phi = matrix.col((1, 0, 0))
        kappa = matrix.col((0.914, 0.279, -0.297))

        Rphi = phi.axis_and_angle_as_r3_rotation_matrix(phi_value, deg=True)
        Rkappa = kappa.axis_and_angle_as_r3_rotation_matrix(kappa_value, deg=True)
        fixed = Rkappa * Rphi

        return self._goniometer_factory.make_goniometer(axis, fixed)

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
        print(FormatCBFMiniPilatusDLS6MSN126.understand(arg))
