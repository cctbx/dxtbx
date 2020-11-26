from __future__ import absolute_import, division, print_function

import sys

from dxtbx import IncorrectFormatError, _ensure_h5_plugin_path
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage


class FormatHDF5(FormatMultiImage, Format):
    def __init__(self, image_file, **kwargs):
        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        # Before opening anything HDF5, ensure that the plugin path is configured
        _ensure_h5_plugin_path()
        FormatMultiImage.__init__(self, **kwargs)
        Format.__init__(self, image_file, **kwargs)

    @staticmethod
    def understand(image_file):
        try:
            with FormatHDF5.open_file(image_file, "rb") as fh:
                return fh.read(8) == b"\211HDF\r\n\032\n"
        except IOError:
            return False


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatHDF5.understand(arg))
