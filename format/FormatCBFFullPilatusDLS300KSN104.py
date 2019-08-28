#!/usr/bin/env python
# FormatCBFFullPilatusDLS300KSN104.py
#
#   Copyright (C) 2017 Diamond Light Source, Richard Gildea
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import math

from dxtbx.format.FormatCBFFullPilatus import FormatCBFFullPilatus
from dxtbx.masking import GoniometerMaskerFactory


class FormatCBFFullPilatusDLS300KSN104(FormatCBFFullPilatus):
    """An image reading class for full CBF format images from Pilatus
    detectors. For DLS I19-2"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        # this depends on DIALS for the goniometer shadow model; if missing
        # simply return False

        header = FormatCBFFullPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS 300K" in record
                and "S/N 3-0104, Diamond" in record
            ):
                return True

        return False

    @staticmethod
    def has_dynamic_shadowing(**kwargs):
        import libtbx

        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if dynamic_shadowing in (libtbx.Auto, "Auto"):
            return False
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)
        super(FormatCBFFullPilatusDLS300KSN104, self).__init__(image_file, **kwargs)

    def get_mask(self, goniometer=None):
        mask = super(FormatCBFFullPilatusDLS300KSN104, self).get_mask()
        if self._dynamic_shadowing:
            gonio_masker = self.get_goniometer_shadow_masker(goniometer=goniometer)
            scan = self.get_scan()
            detector = self.get_detector()
            shadow_mask = gonio_masker.get_mask(detector, scan.get_oscillation()[0])
            assert len(mask) == len(shadow_mask)
            for m, sm in zip(mask, shadow_mask):
                if sm is not None:
                    m &= ~sm
        return mask

    def get_goniometer_shadow_masker(self, goniometer=None):
        if goniometer is None:
            goniometer = self.get_goniometer()

        return GoniometerMaskerFactory.diamond_anvil_cell(
            goniometer, cone_opening_angle=2 * 38 * math.pi / 180
        )


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatCBFFullPilatusDLS300KSN104.understand(arg))
