from __future__ import absolute_import, division, print_function

import sys

from dxtbx.sequence import sequence_factory


def print_sequence(list_of_images):

    s = sequence_factory.sequence(list_of_images, check_headers=True)

    print(s.get_detector())
    print(s.get_beam())
    print(s.get_goniometer())
    print(s.get_scan())


if __name__ == "__main__":
    print_sequence(sys.argv[1:])
