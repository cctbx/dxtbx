"""Print out the resolution (two-theta) of the corners of the detector"""


import math
import sys

from scitbx import matrix

import dxtbx


def resolution_corners(frame):
    """Compute the resolution limit corresponding to the corners of the detector
    surface."""

    detector = frame.get_detector()
    beam = frame.get_beam()

    nfast, nslow = map(int, detector.get_image_size())
    dfast, dslow = detector.get_pixel_size()
    F = matrix.col(detector.get_fast_axis())
    S = matrix.col(detector.get_slow_axis())
    origin = matrix.col(detector.get_origin())

    s0 = -1 * matrix.col(beam.get_sample_to_source_direction())

    for ds in 0, 1:
        for df in 0, 1:
            corner = origin + nfast * dfast * F * df + nslow * dslow * S * ds
            theta = 0.5 * corner.angle(s0)
            print("%.3f" % (beam.get_wavelength() / (2 * math.sin(theta))))


if __name__ == "__main__":
    resolution_corners(dxtbx.load(sys.argv[1]))
