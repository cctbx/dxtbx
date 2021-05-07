"""
An implementation of the CBF image reader for Pilatus images, from the Pilatus
6M SN 114 currently on Diamond VMXi.
"""

import sys

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus


class FormatCBFMiniPilatusDLS6MSN114(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 114 @ DLS."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        year = 0

        for record in header.split("\n"):
            if "# 20" in record:
                year = int(record.replace("-", " ").replace("/", " ").split()[1])
                break

        if year <= 0:
            return False

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 60-0114" in header
            ):
                if year >= 2017:
                    return True
                else:
                    return False

        return False

    def _goniometer(self):
        return self._goniometer_factory.known_axis((0, 1, 0))


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDLS6MSN114.understand(arg))
