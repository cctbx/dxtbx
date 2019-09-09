#!/usr/bin/env python
#   Copyright (C) 2018 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
"""An implementation of the CBF image reader for Eiger images"""

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger


class FormatCBFMiniEigerPetraP14(FormatCBFMiniEiger):
    """A class for reading mini CBF format Eiger images, and correctly
    constructing a model for the experiment from this. This tuned for Petra P14"""

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
        return self._goniometer_factory.known_axis((0, 1, 0))


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatCBFMiniEigerPetraP14.understand(arg))
