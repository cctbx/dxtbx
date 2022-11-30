from __future__ import annotations

import h5py

import dxtbx.nexus
from dxtbx.format.FormatNXmx import FormatNXmx


class FormatNXmxEigerFilewriter(FormatNXmx):

    _cached_file_handle = None

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            if "/entry/instrument/detector/detectorSpecific/eiger_fw_version" in handle:
                return True
        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        super().__init__(image_file, **kwargs)

    def _get_nxmx(self, fh: h5py.File):
        nxmx = dxtbx.nexus.nxmx.NXmx(fh)
        nxentry = nxmx.entries[0]

        nxdetector = nxentry.instruments[0].detectors[0]
        if nxdetector.underload_value is None:
            nxdetector.underload_value = 0

        # data_size is reversed - we should probably be more specific in when
        # we do this, i.e. check data_size is in a list of known reversed
        # values
        for module in nxdetector.modules:
            module.data_size = module.data_size[::-1]
        return nxmx
