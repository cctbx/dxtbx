import h5py

# static enumerations for useless file types

from enum import Enum
HDF5_FileType = Enum("HDF5_FileType", "DATA META MPCCD NXS UNKNOWN")


def get_hdf5_file_type(hdf5_filename):
    """Open the named file, decide if it is probably not a useful one for
    processing..."""

    with h5py.File(hdf5_filename, "r") as f:
        top_level = set(f.keys())
        # case HDF5 data file
        if top_level == {"data"}:
            return HDF5_FileType.DATA

        # case Diamond "meta" file
        if "config" in top_level and "frame_written" in top_level:
            return HDF5_FileType.META

        # NeXus type files contain only /entry at the top (at the moment)
        if top_level == {"entry"}:
            return HDF5_FileType.NXS

        # SACLA MPCCD contains "metadata" and ... other stuff
        if "metadata" in top_level:
            return HDF5_FileType.MPCCD

    return HDF5_FileType.UNKNOWN
