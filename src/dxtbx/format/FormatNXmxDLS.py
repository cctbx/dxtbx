import datetime
import os
from pathlib import Path
from typing import Union

try:
    from functools import cached_property
except ImportError:
    # Python 3.7 compatibility
    # Defined cached_property decorator as a noop
    from dxtbx.nexus.nxmx import cached_property

import h5py

import dxtbx.nexus
from dxtbx.format.FormatNXmx import FormatNXmx


def get_bit_depth_from_meta(meta_file_name):
    with h5py.File(meta_file_name, "r") as f:
        return int(f["/_dectris/bit_depth_image"][()])


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
            if name and name.lower() in ("i03", "i04", "i24", "vmxi"):
                return True
            if name and name.upper().startswith(("DIAMOND BEAMLINE", "DLS")):
                return True
        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        self._legacy = None
        super().__init__(image_file, **kwargs)

    def _start(self):
        super()._start()

        if self._legacy:
            # Get the bit depth from the meta.h5 in order to distinguish masked and
            # saturated pixels. Ideally we would get this from
            # /entry/instrument/detector/bit_depth_readout.
            # See https://jira.diamond.ac.uk/browse/MXGDA-3674
            try:
                self._bit_depth_readout = get_bit_depth_from_meta(self._meta)
            except Exception:
                self._bit_depth_readout = 16

    @cached_property
    def _meta(self):
        return find_meta_filename(self._image_file)

    def _get_nxmx(self, fh: h5py.File):
        nxmx = dxtbx.nexus.nxmx.NXmx(fh)
        nxentry = nxmx.entries[0]

        if self._legacy is None:
            name = dxtbx.nexus.nxmx.h5str(FormatNXmx.get_instrument_name(fh))
            if nxentry.start_time and "I03" in name.upper():
                self._legacy = nxentry.start_time.replace(
                    tzinfo=None
                ) < datetime.datetime(2021, 9, 10, 13, 12, 0)
            elif nxentry.start_time and "I04" in name.upper():
                self._legacy = nxentry.start_time.replace(
                    tzinfo=None
                ) < datetime.datetime(2021, 9, 14, 15, 3, 0)
            elif "VMXI" in name.upper():
                self._legacy = True
            else:
                self._legacy = False

        nxdetector = nxentry.instruments[0].detectors[0]
        if nxdetector.underload_value is None:
            nxdetector.underload_value = 0

        if self._legacy:
            # data_size was reversed, see https://jira.diamond.ac.uk/browse/MXGDA-3676
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
                # At some point the pixel mask was stored with the dimensions reversed -> ignore
                # https://jira.diamond.ac.uk/browse/MXGDA-3675
                nxdetector.pixel_mask = None
        return nxmx
