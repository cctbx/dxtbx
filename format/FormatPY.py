"""Implementation of a base class to read a pickled Python dictionary."""

from __future__ import absolute_import, division, print_function

import sys
import six

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format
import pickle


class FormatPY(Format):
    """Let's take an educated guess as to how to recognize a Python
    pickle file containing a dictionary.  Not easy because there are
    four pickle protocols in Python.  The lowest pickle format only
    gives us two unique bytes by which to recognize a dictionary,
    so also check if we can unpickle the file."""

    @staticmethod
    def understand(image_file):
        try:
            with FormatPY.open_file(image_file, "rb") as fh:
                tag = fh.read(4)
                if (
                    tag[0:2] == b"(d"
                    or tag[0:2] == b"}q"
                    or tag[0:4] == b"\x80\x02}q"
                    or tag[0:4] == b"\x80\x04\x95\xea"
                ):
                    fh.seek(0)
                    if six.PY3:
                        pickle.load(fh, encoding="bytes")
                    else:
                        pickle.load(fh)
                else:
                    return False
        except (IOError, pickle.UnpicklingError):
            return False
        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        super(FormatPY, self).__init__(image_file, **kwargs)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPY.understand(arg))
