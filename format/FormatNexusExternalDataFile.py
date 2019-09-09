from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatHDF5 import FormatHDF5


def find_entries(nx_file):
    """
    Find NXmx entries

    """
    if "entry" in nx_file:
        entry = nx_file["entry"]
        if "NX_class" in entry.attrs:
            if entry.attrs["NX_class"] == "NXentry":
                if "definition" not in entry:
                    return entry
    return None


def is_nexus_external_data_file(filename):
    """
    A hacky function to check if this is a nexus file

    """
    import h5py

    # Get the file handle
    handle = h5py.File(filename, "r")

    # Find the NXmx entries
    entry = find_entries(handle)
    if entry is not None:
        if "instrument" not in entry:
            return True
    return False


class FormatNexusExternalDataFile(FormatHDF5):
    @staticmethod
    def understand(image_file):
        try:
            return is_nexus_external_data_file(image_file)
        except IOError:
            return False

    @classmethod
    def ignore(cls):
        return True


if __name__ == "__main__":
    import sys

    f = FormatNexusExternalDataFile(sys.argv[1])
