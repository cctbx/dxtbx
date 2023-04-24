"""An implementation of the CBF image reader for Eiger images"""


from __future__ import annotations

import sys

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger


class FormatCBFMiniEigerPetraP13(FormatCBFMiniEiger):
    """A class for reading mini CBF format Eiger images, and correctly
    constructing a model for the experiment from this. This tuned for Petra P13"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Eiger mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniEiger.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# detector" in record.lower()
                and "eiger" in record.lower()
                and "E-32-0107" in record
            ):
                return True

        return False

    def _goniometer(self):
        return self._goniometer_factory.known_axis((1, 0, 0))


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniEigerPetraP13.understand(arg))
