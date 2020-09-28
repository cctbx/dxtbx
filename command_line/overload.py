from __future__ import absolute_import, division, print_function

import sys

from dxtbx import load


def overload(image_file):
    i = load(image_file)
    data = i.get_raw_data()
    if not isinstance(data, tuple):
        data = (data,)
    detector = i.get_detector()
    for pid, (d, p) in enumerate(zip(data, detector)):
        if max(d) > p.get_trusted_range()[1]:
            return True
    return False


def run(args=None):
    for image_file in args or sys.argv[1:]:
        if overload(image_file):
            print(image_file)


if __name__ == "__main__":
    run()
