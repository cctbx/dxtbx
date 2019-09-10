from __future__ import absolute_import, division, print_function

import sys

from dxtbx.sweep_filenames import find_matching_images


def print_matching_images(image):
    matching_images = find_matching_images(image)
    for mi in matching_images:
        print(mi)


if __name__ == "__main__":
    print_matching_images(sys.argv[1])
