from __future__ import absolute_import, division, print_function

from dxtbx.format.Format import Format
from libtbx.phil import parse
import os

class FormatPhilLocator(Format):
    @staticmethod
    def understand(image_file):
        if os.path.splitext(image_file)[1] != '.loc': return False
        try:
            return parse(file_name=image_file) is not None
        except Exception:
            return False

if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        print(FormatPhilLocator.understand(arg))
