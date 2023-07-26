"""An implementation of the CBF image reader for Eiger images"""


from __future__ import annotations

import datetime
import sys

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger


class FormatCBFMiniEigerPetraP14(FormatCBFMiniEiger):
    """A class for reading mini CBF format Eiger images, and correctly
    constructing a model for the experiment from this. This tuned for Petra P14"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Eiger mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniEiger.get_cbf_header(image_file)

        # Valid from 22nd May 2021
        expected_serial = "E-32-0129"
        if timestamp := FormatCBFMiniEiger._get_timestamp_from_raw_header(header):
            # We have a timestamp. Let's see what detector we should expect

            # Before 22nd May 2021
            if timestamp < datetime.datetime(2021, 5, 22):
                expected_serial = "E-32-0107"

        # Find the line recording detector serial, and check
        for record in header.split("\n"):
            if (
                "# detector" in record.lower()
                and "eiger" in record.lower()
                and expected_serial in record
            ):
                return True

        return False

    def _goniometer(self):
        return self._goniometer_factory.known_axis((0, 1, 0))


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniEigerPetraP14.understand(arg))
