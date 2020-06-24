from __future__ import absolute_import, division, print_function

import h5py

from dxtbx.format.FormatNexusEigerDLS import FormatNexusEigerDLS


class FormatNexusEigerDLS4M(FormatNexusEigerDLS):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            name = FormatNexusEigerDLS.get_instrument_name(handle)
            if name is None or name.lower() not in (b"vmxi"):
                return False

        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        super(FormatNexusEigerDLS4M, self).__init__(image_file, **kwargs)
        self._dynamic_shadowing = False

    def get_goniometer_shadow_masker(self, goniometer=None):
        return None
