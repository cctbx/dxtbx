import copy
import pickle
import sys

from iotbx.detectors.npy import image_dict_to_unicode

from dxtbx.format.FormatPYunspecified import (
    FormatPYunspecified,
    FormatPYunspecifiedInMemory,
)
from dxtbx.format.FormatStill import FormatStill


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
                data = image_dict_to_unicode(pickle.load(fh, encoding="bytes"))
        except OSError:
            return False

        if "OSC_START" not in data or "OSC_RANGE" not in data:
            return True

        return data["OSC_RANGE"] <= 0


class FormatPYunspecifiedStillInMemory(FormatStill, FormatPYunspecifiedInMemory):
    """Overrides the Format object's init method to accept an image dictionary
    instead of a file name. Used with XFELs when it is desirable to never write
    a file to disk, but to process it only in memory.
    """

    @staticmethod
    def understand(image_file):
        """If it's an image dictionary, we understand this"""
        data = image_file
        if not isinstance(data, dict):
            return False
        if "OSC_START" not in data or "OSC_RANGE" not in data:
            return True
        return data["OSC_RANGE"] <= 0

    def __init__(self, data, **kwargs):
        """@param data In memory image dictionary, alredy initialized"""
        super().__init__(data, **kwargs)
        self._image_file = copy.deepcopy(data)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPYunspecifiedStill.understand(arg))
