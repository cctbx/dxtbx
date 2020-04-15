"""Implementation of a base class to read a pickled Python dictionary."""

from __future__ import absolute_import, division, print_function

import sys

from dxtbx import IncorrectFormatError
from dxtbx.format.FormatFile import FormatFile


class FormatPY(FormatFile):
    """Let's take an educated guess as to how to recognize a Python
    pickle file containing a dictionary.  Not easy because there are
    three pickle protocols in Python 2.7.  Dangerous because the lowest
    pickle format only gives us two unique bytes by which to recognize
    a dictionary.  Could possibly conflict with other image formats."""

    @staticmethod
    def understand(image_file):
        try:
            with FormatPY.open_file(image_file, "rb") as fh:
                tag = fh.read(4)
            return tag[0:2] == b"(d" or tag[0:2] == b"}q" or tag[0:4] == b"\x80\x02}q"
        except IOError:
            return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        super(FormatPY, self).__init__(image_file, **kwargs)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPY.understand(arg))
