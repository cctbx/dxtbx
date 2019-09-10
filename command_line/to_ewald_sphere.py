from __future__ import absolute_import, division, print_function

import sys
from builtins import range

from dxtbx import ImageToEwaldSphere
from dxtbx.sweep import SweepFactory


def to_ewald_sphere(list_of_images):

    sweep = SweepFactory.sweep(list_of_images)
    beam = sweep.get_beam()
    detector = sweep.get_detector()
    gonio = sweep.get_goniometer()
    scan = sweep.get_scan()
    start, end = sweep.get_array_range()
    image_to_ewald_sphere = ImageToEwaldSphere(beam, detector, gonio, scan)

    for frame in range(start, end):
        intensity = sweep[frame]
        x_list = image_to_ewald_sphere(frame)
        for k, x in enumerate(x_list):
            i = intensity[k]
            print("%f %f %f %d" % (x[0], x[1], x[2], i))


if __name__ == "__main__":
    to_ewald_sphere(sys.argv[1:])
