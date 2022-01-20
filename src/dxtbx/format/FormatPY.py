"""Implementation of a base class to read a pickled Python dictionary."""


from __future__ import annotations

import pickle
import sys

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format, abstract


@abstract
class FormatPY(Format):
    """Let's take an educated guess as to how to recognize a Python
    pickle file containing a dictionary.  Not easy because there are
    six pickle protocols in Python.  The lowest pickle format only
    gives us two unique bytes by which to recognize a dictionary,
    so also check if we can unpickle the file."""

    @staticmethod
    def understand(image_file):
        try:
            with FormatPY.open_file(image_file, "rb") as fh:
                tag = fh.read(4)
                if (  # pickle protocols 0-5 where the pickled object is a dictionary
                    tag[0:2] == b"(d"
                    or tag[0:2] == b"}q"
                    or tag[0:4] == b"\x80\x02}q"
                    or tag[0:4] == b"\x80\x03}q"
                    or tag[0:3] == b"\x80\x04\x95"
                    or tag[0:3] == b"\x80\x05\x95"
                ):
                    fh.seek(0)
                    pickle.load(fh, encoding="bytes")
                else:
                    return False
        except (OSError, pickle.UnpicklingError):
            return False
        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        super().__init__(image_file, **kwargs)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPY.understand(arg))
