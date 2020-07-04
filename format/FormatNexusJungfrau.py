from __future__ import absolute_import, division, print_function

import h5py

from dxtbx.format.FormatNexus import FormatNexus


class FormatNexusJungfrau(FormatNexus):
    @staticmethod
    def understand(image_file):
        with h5py.File(image_file, "r") as handle:
            if b"/entry/instrument/JF16M" in handle:
                return True
        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        super(FormatNexusJungfrau, self).__init__(image_file, **kwargs)

    def get_detector(self, index=None):
        detector = self._detector()

        # reassign default values if NULL found in parent class
        # N.B. these came from implementation in nexus.py
        for panel in detector:
            trusted = panel.get_trusted_range()
            if trusted == (-0x7FFFFFFF, 0x7FFFFFFF):
                trusted = -400, 90000
            panel.set_trusted_range(trusted)

        return detector
