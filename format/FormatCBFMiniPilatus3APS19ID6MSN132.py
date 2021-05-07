import sys

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus


class FormatCBFMiniPilatus3AOS19ID6MSN132(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus3 images for 6M SN 132 @ APS19ID."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS3" in record
                and "S/N 60-0132" in header
            ):
                return True

        return False

    def _goniometer(self):
        """19ID has reversed goniometer"""

        return self._goniometer_factory.single_axis_reverse()


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatus3AOS19ID6MSN132.understand(arg))
