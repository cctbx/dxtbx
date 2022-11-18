from __future__ import annotations

import h5py
import numpy as np

import dxtbx.nexus
from dxtbx.format.FormatNXmx import FormatNXmx


class FormatNXmxEiger(FormatNXmx):

    _cached_file_handle = None

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            name = dxtbx.nexus.nxmx.h5str(FormatNXmxEiger.get_instrument_name(handle))
            # known back to front data size examples
            if name and name.lower() in ("biomax",):
                return True
        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        self._legacy = None
        super().__init__(image_file, **kwargs)

    def _start(self):
        super()._start()

    def _get_nxmx(self, fh: h5py.File):
        nxmx = dxtbx.nexus.nxmx.NXmx(fh)
        nxentry = nxmx.entries[0]

        nxdetector = nxentry.instruments[0].detectors[0]
        if nxdetector.underload_value is None:
            nxdetector.underload_value = 0

        if self._legacy is None:
            name = dxtbx.nexus.nxmx.h5str(FormatNXmx.get_instrument_name(fh)).lower()
            if name in ("biomax",):
                self._legacy = np.all(nxdetector.modules[0].data_size == (4150, 4371))

        if self._legacy:
            for module in nxdetector.modules:
                module.data_size = module.data_size[::-1]
            try:
                nxdetector.pixel_mask
            except KeyError:
                nxdetector.pixel_mask = None
            if (
                nxdetector.pixel_mask is not None
                and nxdetector.pixel_mask.shape
                != tuple(nxdetector.modules[0].data_size)
            ):
                nxdetector.pixel_mask = None
        return nxmx
