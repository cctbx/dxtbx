#!/usr/bin/env python
# FormatCBFMiniPilatusDLS6MSN100.py
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images, from the Pilatus
# 6M SN 100 currently on Diamond I04.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.format.FormatStill import FormatStill


class FormatCBFMiniPilatusDLS6MSN114DMM(FormatCBFMiniPilatus, FormatStill):
    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        import os

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

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatCBFMiniPilatus.__init__(self, image_file, **kwargs)
        FormatStill.__init__(self, image_file, **kwargs)

        return

    def _goniometer(self):
        return None

    def _scan(self):
        return None


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDLS6MSN114DMM.understand(arg))
