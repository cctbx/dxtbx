#!/usr/bin/env python
# FormatCBFMiniPilatusSPring8_6MSN125.py
#
#  Copyright (C) (2015) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus

# from dxtbx.model import ParallaxCorrectedPxMmStrategy
# from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask


class FormatCBFMiniPilatusSPring8_6MSN125(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 125, normally
    at Spring8 BL41XU"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 60-0125" in header
            ):
                return True

        return False

    def _goniometer(self):
        """Return a model for a simple single-axis reversed direction goniometer."""

        return self._goniometer_factory.single_axis_reverse()
