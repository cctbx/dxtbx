from __future__ import annotations

import re

import h5py
import nxmx
from packaging import version

from dxtbx.format.FormatNXmxEigerFilewriter import FormatNXmxEigerFilewriter

class FormatNXmxEigerFilewriterCHESS(FormatNXmxEigerFilewriter):
    _cached_file_handle = None

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            if "/entry/instrument/detector/detector_number" in handle:
                if (
                    handle["/entry/instrument/detector/detector_number"][()]
                    == b'E-32-0123'
                ):
                    return True
        return False

    def _goniometer(self):
        return self._goniometer_factory.known_axis((-1, 0, 0))
