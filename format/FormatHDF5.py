from __future__ import absolute_import, division, print_function

import sys

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format
from dxtbx.format.FormatHDF5Exceptions import HDF5_FileType, get_hdf5_file_type
from dxtbx.format.FormatMultiImage import FormatMultiImage


class FormatHDF5(FormatMultiImage, Format):
    def __init__(self, image_file, **kwargs):
        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatMultiImage.__init__(self, **kwargs)
        Format.__init__(self, image_file, **kwargs)

    @staticmethod
    def understand(image_file):
        try:
            filetype = get_hdf5_file_type(image_file)
            return filetype in (HDF5_FileType.NXS, HDF5_FileType.MPCCD)
        except (IOError, OSError):
            return False


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatHDF5.understand(arg))
