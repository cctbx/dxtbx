from __future__ import absolute_import, division, print_function

import sys

from dxtbx.sweep import sweep_factory


def print_sweep(list_of_images):

    s = sweep_factory.sweep(list_of_images, check_headers=True)

    print(s.get_detector())
    print(s.get_beam())
    print(s.get_goniometer())
    print(s.get_scan())


if __name__ == "__main__":
    print_sweep(sys.argv[1:])
