from __future__ import absolute_import, division, print_function

from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage


class FormatHDF5(FormatMultiImage, Format):
    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
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
    import sys

    for arg in sys.argv[1:]:
        print(FormatHDF5.understand(arg))
