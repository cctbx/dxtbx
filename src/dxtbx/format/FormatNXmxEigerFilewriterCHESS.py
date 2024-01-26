from __future__ import annotations

import re

import h5py
import nxmx
from packaging import version

from dxtbx.format.FormatNXmxEigerFilewriter import FormatNXmxEigerFilewriter

DATA_FILE_RE = re.compile(r"data_\d{6}")


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

    def _get_nxmx(self, fh: h5py.File):
        nxmx_obj = nxmx.NXmx(fh)
        nxentry = nxmx_obj.entries[0]

        nxdetector = nxentry.instruments[0].detectors[0]
        if nxdetector.underload_value is None:
            nxdetector.underload_value = 0

        # older firmware versions had the detector dimensions inverted
        fw_version_string = (
            fh["/entry/instrument/detector/detectorSpecific/eiger_fw_version"][()]
            .decode()
            .replace("release-", "")
        )
        
        #print(f'detected Eiger firmware version: {fw_version_string}')
        if version.parse(fw_version_string) < version.parse("2022.1.2"):
            #print('swapping module size')
            for module in nxdetector.modules:
                module.data_size = module.data_size[::-1]
        return nxmx_obj
    
    def _goniometer(self):
        return self._goniometer_factory.known_axis((-1, 0, 0))