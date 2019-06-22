from __future__ import absolute_import, division, print_function

from builtins import range
import math
import os
import random

import dxtbx
import pytest
import six.moves.cPickle as pickle
from cctbx.eltbx import attenuation_coefficient
from dxtbx.model import ParallaxCorrectedPxMmStrategy
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.array_family import flex


@pytest.fixture
def model(dials_regression):
    filename = os.path.join(dials_regression, "image_examples", "XDS", "XPARM.XDS")

    models = dxtbx.load(filename)
    detector = models.get_detector()
    assert len(detector) == 1
    detector = detector[0]
    t0 = 0.320

    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(models.get_beam().get_wavelength()) / 10.0
    pixel_size = detector.get_pixel_size()

    return {"mu": mu, "t0": t0, "detector": detector, "pixel_size": pixel_size}


def correct_gold(model, xy):
    mu = model["mu"]
    t0 = model["t0"]
    s1 = matrix.col(model["detector"].get_lab_coord(xy)).normalize()
    d0 = matrix.col(model["detector"].get_origin())
    d1 = matrix.col(model["detector"].get_fast_axis())
    d2 = matrix.col(model["detector"].get_slow_axis())
    dn = d1.cross(d2)
    cos_theta = s1.dot(dn)
    t = t0 / cos_theta
    o = (1.0 / mu) - (t + 1.0 / mu) * math.exp(-mu * t)
    cx = xy[0] + s1.dot(d1) * o
    cy = xy[1] + s1.dot(d2) * o
    return (cx / model["pixel_size"][0], cy / model["pixel_size"][1])


def test_correction_on_random_coordinates(model):
    convert = ParallaxCorrectedPxMmStrategy(model["mu"], model["t0"])
    for i in range(10000):
        xy = matrix.col((random.uniform(-1000, 1000), random.uniform(-1000, 1000)))

        # Do the forward and reverse corrections
        xy_corr_gold = matrix.col(correct_gold(model, xy))
        xy_corr = matrix.col(convert.to_pixel(model["detector"], xy))
        xy_corr_inv = matrix.col(convert.to_millimeter(model["detector"], xy_corr))

        # Check the values
        assert abs(xy_corr_gold - xy_corr) < 1e-7
        assert abs(xy_corr_inv - xy) < 1e-3


def test_array(model):
    random_coord = lambda: (random.uniform(-1000, 1000), random.uniform(-1000, 1000))
    xy = flex.vec2_double([random_coord() for i in range(100)])
    xy_corr = model["detector"].get_lab_coord(xy)
    xy_corr_panel = model["detector"].get_lab_coord(xy)
    xy_corr_gold = [model["detector"].get_lab_coord(xy_single) for xy_single in xy]
    assert approx_equal(xy_corr, xy_corr_gold)
    assert approx_equal(xy_corr_panel, xy_corr_gold)


def test_inverted_axis():
    def get_values(invert_y):
        from dxtbx.model.beam import BeamFactory

        beam = BeamFactory.simple(wavelength=1)

        if invert_y:
            y_direction = "-y"
        else:
            y_direction = "+y"

        from dxtbx.model.detector import DetectorFactory

        detector = DetectorFactory.simple(
            sensor=DetectorFactory.sensor("PAD"),
            distance=100,
            beam_centre=[50, 50],
            fast_direction="+x",
            slow_direction=y_direction,
            pixel_size=[0.1, 0.1],
            image_size=[1000, 1000],
        )[0]

        wavelength = beam.get_wavelength()
        thickness = 0.5
        table = attenuation_coefficient.get_table("Si")
        mu = table.mu_at_angstrom(wavelength) / 10.0
        t0 = thickness

        for panel in detector:
            panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
        v1 = detector.pixel_to_millimeter((0, 0))
        v2 = detector.pixel_to_millimeter((1000, 1000))

        return v1, v2

    v11, v12 = get_values(False)
    v21, v22 = get_values(False)

    assert abs(v11[0] - v21[0]) < 1e-7
    assert abs(v11[1] - v21[1]) < 1e-7
    assert abs(v12[0] - v22[0]) < 1e-7
    assert abs(v12[1] - v22[1]) < 1e-7


def test_offset_px_mm_strategy():
    from dxtbx.model import Panel
    from dxtbx.model import OffsetParallaxCorrectedPxMmStrategy

    # for future reference this is the array the same shape
    # as the image in pixels with offsets in pixels

    dx = flex.double(flex.grid(10, 10), 1)
    dy = flex.double(flex.grid(10, 10), 1)

    strategy = OffsetParallaxCorrectedPxMmStrategy(1, 1, dx, dy)

    p = Panel()
    p.set_image_size((10, 10))
    p.set_mu(1)
    p.set_thickness(1)
    p.set_px_mm_strategy(strategy)

    d = p.to_dict()

    pnew = Panel.from_dict(d)
    assert pnew == p

    pnew = pickle.loads(pickle.dumps(pnew))
    assert pnew == p
