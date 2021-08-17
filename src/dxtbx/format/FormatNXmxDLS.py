import os
from pathlib import Path
from typing import Union

import h5py

from dials.array_family import flex

import dxtbx.nexus
from dxtbx.format.FormatNXmx import FormatNXmx

# Hack to switch off the FormatNexusEigerDLS format class
# from dxtbx.format.FormatNexusEigerDLS import FormatNexusEigerDLS
# FormatNexusEigerDLS.understand = lambda image_file: False


def get_bit_depth_from_meta(meta_file_name):
    with h5py.File(meta_file_name, "r") as f:
        return int(f["/_dectris/bit_depth_image"][()][0])


def find_meta_filename(master_like: Union[str, Path]) -> Union[str, Path]:
    """
    Find the path to the '..._meta.h5' file in the same directory as the master file.

    Args:
        master_like:  File path of the master HDF5 file.

    Returns:
        File path of the HDF5 metadata '..._meta.h5' file.
    """

    def _local_visit(name):
        obj = f[name]
        if not hasattr(obj, "keys"):
            return None
        for k in obj.keys():
            kclass = obj.get(k, getlink=True, getclass=True)
            if kclass is h5py._hl.group.ExternalLink:
                kfile = obj.get(k, getlink=True).filename
                if kfile.split(".")[0].endswith("meta"):
                    return kfile

    master_dir = os.path.split(master_like)[0]
    with h5py.File(master_like) as f:
        meta_filename = f.visit(_local_visit)

    return os.path.join(master_dir, meta_filename)


class FormatNXmxDLS(FormatNXmx):

    _cached_file_handle = None

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file, "r") as handle:
            name = dxtbx.nexus.nxmx.h5str(FormatNXmxDLS.get_instrument_name(handle))
            if name is None:
                return False
            if name.lower() in {"i03", "i04", "i24", "vmxi"}:
                return True
            if name.upper().startswith("DLS ") and "i19-2" not in name.lower():
                return True
        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        super().__init__(image_file, **kwargs)
        # Get the bit depth from the meta.h5 in order to distinguish masked and
        # saturated pixels. Ideally we would get this from
        # /entry/instrument/detector/bit_depth_readout.
        # See https://jira.diamond.ac.uk/browse/MXGDA-3674
        try:
            meta = find_meta_filename(image_file)
            self._bit_depth_image = get_bit_depth_from_meta(meta)
        except Exception:
            self._bit_depth_image = 16

    def _start(self):
        super()._start()
        # Due to a bug the dimensions (but not the values) of the pixel_mask array
        # are reversed. See https://jira.diamond.ac.uk/browse/MXGDA-3675.
        if self._static_mask:
            for m in self._static_mask:
                m.reshape(flex.grid(reversed(m.all())))
        # /entry/instrument/detector/module/data_size is also reversed:
        # https://jira.diamond.ac.uk/browse/MXGDA-3676
        for panel in self._detector_model:
            panel.set_image_size(tuple(reversed(panel.get_image_size())))

    def get_raw_data(self, index):
        if self._cached_file_handle is None:
            self._cached_file_handle = h5py.File(self._image_file, "r")

        # /entry/instrument/detector/module/data_size is reversed:
        # https://jira.diamond.ac.uk/browse/MXGDA-3676
        nxmx = dxtbx.nexus.nxmx.NXmx(self._cached_file_handle)
        nxdata = nxmx.entries[0].data[0]
        nxdetector = nxmx.entries[0].instruments[0].detectors[0]
        for module in nxdetector.modules:
            module.data_size = module.data_size[::-1]
        data = dxtbx.nexus.get_raw_data(nxdata, nxdetector, index)[0]

        if self._bit_depth_image:
            # if 32 bit then it is a signed int, I think if 8, 16 then it is
            # unsigned with the highest two values assigned as masking values
            if self._bit_depth_image == 32:
                top = 2 ** 31
            else:
                top = 2 ** self._bit_depth_image
            d1d = data.as_1d()
            d1d.set_selected(d1d == top - 1, -1)
            d1d.set_selected(d1d == top - 2, -2)
        return data
