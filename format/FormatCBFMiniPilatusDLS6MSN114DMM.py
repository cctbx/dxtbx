"""
An implementation of the CBF image reader for Pilatus images, from the Pilatus
6M SN 100 currently on Diamond I04.
"""


import os
import sys

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.format.FormatStill import FormatStill


class FormatCBFMiniPilatusDLS6MSN114DMM(FormatCBFMiniPilatus, FormatStill):
    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        if os.environ.get("I02_DMM", "0") != "1":
            return False

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 60-0114-F" in header
            ):
                return True

        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        FormatCBFMiniPilatus.__init__(self, image_file, **kwargs)
        FormatStill.__init__(self, image_file, **kwargs)

        return

    def _goniometer(self):
        return None

    def _scan(self):
        return None


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDLS6MSN114DMM.understand(arg))
