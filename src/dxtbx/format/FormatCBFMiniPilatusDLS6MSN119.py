"""
An implementation of the CBF image reader for Pilatus images, from the Pilatus
6M SN 119 currently on Diamond I24.
"""


import sys

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus


class FormatCBFMiniPilatusDLS6MSN119(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 119 @ DLS."""

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
                and "S/N 60-0119" in header
            ):
                if year >= 2015:
                    return True
                else:
                    return False

        return False

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header, though for miniCBF
        there are limited options for this."""

        # FIXME this is currently a mess and needs properly sorting :(

        if "+FAST" in self._cif_header_dictionary.get("Oscillation_axis"):
            # plate mode - though should be -FAST
            return self._goniometer_factory.known_axis((-1, 0, 0))
        elif "-FAST" in self._cif_header_dictionary.get("Oscillation_axis"):
            # plate mode - after the fix
            return self._goniometer_factory.known_axis((-1, 0, 0))
        return self._goniometer_factory.known_axis((0, 1, 0))


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDLS6MSN119.understand(arg))
