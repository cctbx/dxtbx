from __future__ import absolute_import, division, print_function

from future import standard_library
from dxtbx import IncorrectFormatError
import sys

standard_library.install_aliases()

import copy
import pickle
import six
from dxtbx.format.FormatStill import FormatStill
from dxtbx.format.FormatPYunspecified import FormatPYunspecified
from dxtbx.format.FormatPYunspecified import FormatPYunspecifiedInMemory


class FormatPYunspecifiedStill(FormatStill, FormatPYunspecified):
    @staticmethod
    def understand(image_file):
        """Seems like the static method wastes a lot of effort here; it's not possible to
        just read the first few bytes; instead understand() reads the entire first data
        item in the file; an entire binary image.  This data is then read again in the
        _start() method and again in the detectorbase constructor.
        """

        try:
            with FormatPYunspecified.open_file(image_file, "rb") as fh:
                if six.PY3:
                    data = pickle.load(fh, encoding="bytes")  # lgtm
                    # the '# lgtm' comment disables an lgtm false positive and
                    # can be removed once we move to Python 3 only
                    data = {key.decode("ascii"): value for key, value in data.items()}
                else:
                    data = pickle.load(fh)
        except IOError:
            return False

        if "OSC_START" not in data or "OSC_RANGE" not in data:
            return True

        return data["OSC_RANGE"] <= 0

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatPYunspecified.__init__(self, image_file, **kwargs)


class FormatPYunspecifiedStillInMemory(FormatStill, FormatPYunspecifiedInMemory):
    """Overrides the Format object's init method to accept an image dictionary
    instead of a file name. Used with XFELs when it is desirable to never write
    a file to disk, but to process it only in memory.
    """

    @staticmethod
    def understand(image_file):
        """ If it's an image dictionary, we understand this """
        data = image_file
        if not isinstance(data, dict):
            return False
        if "OSC_START" not in data or "OSC_RANGE" not in data:
            return True
        return data["OSC_RANGE"] <= 0

    def __init__(self, data, **kwargs):
        """ @param data In memory image dictionary, alredy initialized """
        FormatPYunspecifiedInMemory.__init__(self, data, **kwargs)

        self._image_file = copy.deepcopy(data)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPYunspecifiedStill.understand(arg))
