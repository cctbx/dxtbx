from __future__ import annotations

import datetime
from functools import cached_property
from pathlib import Path

import h5py
import numpy as np
import nxmx

from dxtbx.format.FormatNXmx import FormatNXmx


def get_bit_depth_from_meta(meta_file_name):
    with h5py.File(meta_file_name) as f:
        return int(f["/_dectris/bit_depth_image"][()])


def find_meta_filename(master_like: Path) -> Path:
    """
    Find the path to the '..._meta.h5' file in the same directory as the master file.
    Args:
        master_like:  File path of the master HDF5 file.

    Returns:
        File path of the HDF5 metadata '..._meta.h5' file.

    Raises:
        RuntimeError: If no meta-file was found.
    """

    def _local_visit(name) -> Path | None:
        obj = f[name]
        if not hasattr(obj, "keys"):
            return None
        for k in obj.keys():
            kclass = obj.get(k, getlink=True, getclass=True)
            if kclass is h5py._hl.group.ExternalLink:
                kfile = master_like.parent / obj.get(k, getlink=True).filename
                if kfile.stem.endswith("meta"):
                    return kfile
        return None

    master_like = Path(master_like)
    with h5py.File(master_like) as f:
        meta_filename = f.visit(_local_visit)

    if meta_filename is None:
        # Try a fallback of looking for a file named the same but "meta.h5"
        look_for = master_like.with_name(
            master_like.stem.removesuffix("_master") + "_meta.h5"
        )
        if look_for.is_file():
            meta_filename = look_for
        else:
            raise RuntimeError(f"Could not find h5 meta file for {master_like}")

    return meta_filename


class FormatNXmxDLS(FormatNXmx):

    _cached_file_handle = None

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            name = nxmx.h5str(FormatNXmxDLS.get_instrument_name(handle))
            if name and name.lower() in ("i03", "i04", "i24", "vmxi"):
                return True
            if name and "ebic" in name.lower():
                return False
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
            except RuntimeError:
                self._bit_depth_readout = 16

    @cached_property
    def _meta(self) -> Path:
        return find_meta_filename(self._image_file)

    def _get_nxmx(self, fh: h5py.File):
        nxmx_obj = nxmx.NXmx(fh)
        nxentry = nxmx_obj.entries[0]

        nxdetector = nxentry.instruments[0].detectors[0]
        if nxdetector.underload_value is None:
            nxdetector.underload_value = 0

        if self._legacy is None:
            name = nxmx.h5str(FormatNXmx.get_instrument_name(fh))
            if nxentry.start_time and "I03" in name.upper():
                self._legacy = nxentry.start_time.replace(
                    tzinfo=None
                ) < datetime.datetime(2021, 9, 10, 13, 12, 0)
            elif nxentry.start_time and "I04" in name.upper():
                self._legacy = nxentry.start_time.replace(
                    tzinfo=None
                ) < datetime.datetime(2021, 9, 14, 15, 3, 0)
            elif nxentry.start_time and "VMXM" in name.upper():
                # don't have a time stamp for precisely when this happened to VMXm
                self._legacy = nxentry.start_time.replace(
                    tzinfo=None
                ) < datetime.datetime(2022, 1, 1, 1, 0, 0)
            elif "VMXI" in name.upper():
                # Until some point between July 2021 and  November 2021 the
                # data_size recorded in the VMXi NeXus files were reversed
                self._legacy = np.all(nxdetector.modules[0].data_size == (2068, 2162))

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
        return nxmx_obj
