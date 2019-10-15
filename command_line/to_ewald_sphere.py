from __future__ import absolute_import, division, print_function

import sys
from builtins import range

from dxtbx import ImageToEwaldSphere
from dxtbx.sequence import SequenceFactory


def to_ewald_sphere(list_of_images):

    sequence = SequenceFactory.sequence(list_of_images)
    beam = sequence.get_beam()
    detector = sequence.get_detector()
    gonio = sequence.get_goniometer()
    scan = sequence.get_scan()
    start, end = sequence.get_array_range()
    image_to_ewald_sphere = ImageToEwaldSphere(beam, detector, gonio, scan)

    for frame in range(start, end):
        intensity = sequence[frame]
        x_list = image_to_ewald_sphere(frame)
        for k, x in enumerate(x_list):
            i = intensity[k]
            print("%f %f %f %d" % (x[0], x[1], x[2], i))


if __name__ == "__main__":
    to_ewald_sphere(sys.argv[1:])
