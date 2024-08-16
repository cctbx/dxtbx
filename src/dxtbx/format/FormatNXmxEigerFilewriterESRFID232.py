from __future__ import annotations

import re

import h5py

from dxtbx.format.FormatNXmxEigerFilewriter import FormatNXmxEigerFilewriter

DATA_FILE_RE = re.compile(r"data_\d{6}")


class FormatNXmxEigerFilewriterESRFID232(FormatNXmxEigerFilewriter):
    _cached_file_handle = None

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            if "/entry/instrument/detector/detector_number" in handle:
                if (
                    handle["/entry/instrument/detector/detector_number"][()]
                    == b"E-18-0133"
                ):
                    return True
        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        super().__init__(image_file, **kwargs)

    def _start(self):
        super()._start()

    def _goniometer(self):
        return self._goniometer_factory.known_axis((0, 1, 0))
