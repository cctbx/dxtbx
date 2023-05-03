from __future__ import annotations

import re

import h5py
import nxmx

from scitbx.array_family import flex

from dxtbx.format.FormatNXmx import FormatNXmx
from dxtbx.nexus import _dataset_as_flex, get_detector_module_slices

DATA_FILE_RE = re.compile(r"data_\d{6}")


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
        nxmx_obj = nxmx.NXmx(fh)
        nxentry = nxmx_obj.entries[0]

        # possibly useful datum for later
        try:
            self._bit_depth_image = fh["entry/instrument/detector/bit_depth_image"][()]
        except KeyError:
            self._bit_depth_image = None

        nxdetector = nxentry.instruments[0].detectors[0]
        if nxdetector.underload_value is None:
            nxdetector.underload_value = 0

        # data_size is reversed - we should probably be more specific in when
        # we do this, i.e. check data_size is in a list of known reversed
        # values
        for module in nxdetector.modules:
            module.data_size = module.data_size[::-1]
        return nxmx_obj

    def get_raw_data(self, index):
        nxmx_obj = self._get_nxmx(self._cached_file_handle)
        nxdata = nxmx_obj.entries[0].data[0]
        nxdetector = nxmx_obj.entries[0].instruments[0].detectors[0]

        bit_depth = self._bit_depth_image

        raw_data = get_raw_data(nxdata, nxdetector, index, bit_depth)

        if self._bit_depth_readout:
            # if 32 bit then it is a signed int, I think if 8, 16 then it is
            # unsigned with the highest two values assigned as masking values
            if self._bit_depth_readout == 32:
                top = 2**31
            else:
                top = 2**self._bit_depth_readout
            for data in raw_data:
                d1d = data.as_1d()
                d1d.set_selected(d1d == top - 1, -1)
                d1d.set_selected(d1d == top - 2, -2)
        return raw_data


def get_raw_data(
    nxdata: nxmx.NXdata, nxdetector: nxmx.NXdetector, index: int, bit_depth: int
) -> tuple[flex.float | flex.double | flex.int, ...]:
    """Return the raw data for an NXdetector.

    This will be a tuple of flex.float, flex.double or flex.int arrays, of length equal
    to the number of modules. The result is intended to be compatible with the
    get_raw_data() method of dxtbx format classes.

    Specialized version for files written by the DECTRIS EIGER filewriter, which don't
    write a virtual dataset, and instead we have to manually identify the correct
    data_nnnnnn sub-dataset for a given index.
    """
    data_subsets = [v for k, v in sorted(nxdata.items()) if DATA_FILE_RE.match(k)]
    for data in data_subsets:
        if index < data.shape[0]:
            break
        index -= data.shape[0]
    if index >= data.shape[0]:
        raise IndexError(f"Out of range index for raw data {index}")
    all_data = []
    sliced_outer = data[index]
    for module_slices in get_detector_module_slices(nxdetector):
        data_as_flex = _dataset_as_flex(sliced_outer, tuple(module_slices), bit_depth)
        all_data.append(data_as_flex)
    return tuple(all_data)
