import os
import random

from scitbx import matrix

import dxtbx
from dxtbx.model import parallax_correction, parallax_correction_inv


def correct_gold(detector, attlen, xy):
    s1 = matrix.col(detector.get_lab_coord(xy))
    lab = attlen * s1 / s1.length()
    d = detector.get_d_matrix()
    d0 = matrix.col((d[0], d[3], d[6]))
    d1 = matrix.col((d[1], d[4], d[7]))
    mm0 = d0.dot(lab) / d0.length()
    mm1 = d1.dot(lab) / d1.length()
    mmcal = matrix.col((xy[0] + mm0, xy[1] + mm1))
    return mmcal


def test_run(dials_regression):
    filename = os.path.join(dials_regression, "image_examples", "XDS", "XPARM.XDS")

    models = dxtbx.load(filename)
    detector = models.get_detector()
    assert len(detector) == 1
    detector = detector[0]
    attlen = 0.252500934883
    distance = detector.get_distance()
    origin = detector.get_ray_intersection(detector.get_normal())

    for i in range(10000):
        # Generate some random coordinates
        xy = matrix.col((random.uniform(-1000, 1000), random.uniform(-1000, 1000)))

        # Do the forward and reverse corrections
        corr_gold = matrix.col(correct_gold(detector, attlen, xy))
        corr = matrix.col(parallax_correction(distance, attlen, origin, xy))
        corr_inv = matrix.col(parallax_correction_inv(distance, attlen, origin, corr))

        # Check the values
        assert abs(corr_gold - corr) < 1e-7
        assert abs(corr_inv - xy) < 1e-3
