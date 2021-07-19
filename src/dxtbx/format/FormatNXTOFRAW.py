from __future__ import absolute_import, division, print_function

from sys import argv

import h5py

from dxtbx import IncorrectFormatError
from dxtbx.format.FormatHDF5 import FormatHDF5


class FormatNXTOFRAW(FormatHDF5):

    """
    Class to read NXTOFRAW files as defined in
    https://www.nexusformat.org/TOFRaw.html

    See also
    https://manual.nexusformat.org/classes/applications/NXtofraw.html

    """

    def __init__(self, image_file):
        if not FormatNXTOFRAW.understand(image_file):
            raise IncorrectFormatError(self, image_file)

    @staticmethod
    def understand(image_file):
        try:
            return FormatNXTOFRAW.is_nxtofraw_file(image_file)
        except IOError:
            return False

    @staticmethod
    def field_in_file(field, nxs_file):
        def field_in_file_recursive(field, nxs_obj):
            if field in nxs_obj.name.split("/"):
                return True
            else:
                if isinstance(nxs_obj, h5py.Group):
                    for i in nxs_obj.values():
                        if field_in_file_recursive(field, i):
                            return True
                return False

        for i in nxs_file.values():
            if field_in_file_recursive(field, i):
                return True
        return False

    @staticmethod
    def is_nxtofraw_file(image_file):

        """
        Confirms if image_file conforms to NXTOFRAW format
        by checking if required_fields are present

        """

        required_fields = ["detector_1"]

        def required_fields_present(required_fields, image_file):
            with h5py.File(image_file, "r") as handle:
                for i in required_fields:
                    if not FormatNXTOFRAW.field_in_file(i, handle):
                        return False
                return True

        if not FormatHDF5.understand(image_file):
            return False

        return required_fields_present(required_fields, image_file)


if __name__ == "__main__":
    for arg in argv[1:]:
        print(FormatNXTOFRAW.understand(arg))
