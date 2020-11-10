import h5py

# static enumerations for useless file types

HDF5_UNKNOWN = -1
HDF5_NXS_FILE = 0
HDF5_DATA_FILE = 1
HDF5_META_FILE = 2


def hdf5_file_type(hdf5_filename):
    """Open the named file, decide if it is probably not a useful one for
    processing..."""

    with h5py.File(hdf5_filename, "r") as f:
        top_level = set(f.keys())
        # case HDF5 data file
        if top_level == set(["data"]):
            return HDF5_DATA_FILE

        # case Diamond "meta" file
        if "config" in top_level and "frame_written" in top_level:
            return HDF5_META_FILE

        if top_level == set(["entry"]):
            return HDF5_NXS_FILE

    return HDF5_UNKNOWN
