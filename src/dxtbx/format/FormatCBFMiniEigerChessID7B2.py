from __future__ import annotations

import sys

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger


class FormatCBFMiniEigerChessID7B2(FormatCBFMiniEiger):
    """A class for reading mini CBF format Eiger16M images for S/N E-32-0123
    installed at CHESS ID7B2, which has an inverted goniometer axis."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Eiger mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniEiger.get_cbf_header(image_file)
        for record in header.split("\n"):
            if "# Detector: Dectris EIGER2 Si 16M, S/N E-32-0123" in record:
                return True

        return False

    def _goniometer(self):
        return self._goniometer_factory.known_axis((-1, 0, 0))


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniEigerChessID7B2.understand(arg))
