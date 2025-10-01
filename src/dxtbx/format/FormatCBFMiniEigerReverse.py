"""
An implementation of the CBF image reader for Eiger images for beamlines
with a horizontal, but reversed rotation axis.
"""

from __future__ import annotations

import sys

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger


class FormatCBFMiniEigerReverse(FormatCBFMiniEiger):
    """A class for reading mini CBF format Eiger images for beamlines with
    a horizontal reversed rotation axis."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Eiger mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniEiger.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "# Detector" in record and "Eiger" in record:
                if "S/N E-32-0111" in record:
                    # S/N E-32-0111: SSRF BL17U1
                    return True

        return False

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header, though for miniCBF
        there are limited options for this."""

        return self._goniometer_factory.single_axis_reverse()


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniEigerReverse.understand(arg))
