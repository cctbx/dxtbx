import math
import os

import pytest

import libtbx.load_env
from libtbx import easy_pickle
from libtbx.phil import parse
from libtbx.test_utils import Exception_expected, approx_equal
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.model.goniometer import (
    Goniometer,
    GoniometerFactory,
    KappaGoniometer,
    MultiAxisGoniometer,
    goniometer_phil_scope,
)


def _compare_tuples(a, b, tolerance=1.0e-6):
    assert len(a) == len(b)

    for va, vb in zip(a, b):
        assert math.fabs(va - vb) <= tolerance


def test_goniometer():
    """A test class for the goniometer class."""

    axis = (1, 0, 0)
    fixed = (1, 0, 0, 0, 1, 0, 0, 0, 1)

    xg = Goniometer(axis, fixed)

    assert len(xg.get_rotation_axis()) == 3
    assert len(xg.get_fixed_rotation()) == 9

    _compare_tuples(xg.get_rotation_axis(), axis)
    _compare_tuples(xg.get_fixed_rotation(), fixed)

    single = GoniometerFactory.single_axis()

    assert len(single.get_rotation_axis()) == 3
    assert len(single.get_fixed_rotation()) == 9

    _compare_tuples(single.get_rotation_axis(), axis)
    _compare_tuples(single.get_fixed_rotation(), fixed)

    kappa = GoniometerFactory.kappa(50.0, 0.0, 0.0, 0.0, "-y", "omega")

    assert len(kappa.get_rotation_axis()) == 3
    assert len(kappa.get_fixed_rotation()) == 9
    _compare_tuples(kappa.get_rotation_axis(), axis)
    _compare_tuples(kappa.get_fixed_rotation(), fixed)

    kappa = GoniometerFactory.kappa(50.0, 0.0, 0.0, 0.0, "-y", "omega")

    assert len(kappa.get_rotation_axis()) == 3
    assert len(kappa.get_fixed_rotation()) == 9

    _compare_tuples(kappa.get_rotation_axis(), axis)
    _compare_tuples(kappa.get_fixed_rotation(), fixed)

    kappa = GoniometerFactory.kappa(50.0, 0.0, 0.0, 0.0, "-y", "phi")

    assert len(kappa.get_rotation_axis()) == 3
    assert len(kappa.get_fixed_rotation()) == 9

    _compare_tuples(kappa.get_rotation_axis(), axis)
    _compare_tuples(kappa.get_fixed_rotation(), fixed)

    kappa = GoniometerFactory.kappa(50.0, 0.0, 30.0, 0.0, "-y", "omega")

    assert len(kappa.get_rotation_axis()) == 3
    assert len(kappa.get_fixed_rotation()) == 9

    _compare_tuples(kappa.get_rotation_axis(), axis)
    with pytest.raises(AssertionError):
        _compare_tuples(kappa.get_fixed_rotation(), fixed)

    dxtbx_dir = libtbx.env.dist_path("dxtbx")

    image = os.path.join(dxtbx_dir, "tests", "phi_scan_001.cbf")
    assert GoniometerFactory.imgCIF(image)

    kappa = GoniometerFactory.kappa(50.0, -10.0, 30.0, 0.0, "-y", "phi")

    s = easy_pickle.dumps(kappa)
    kappa2 = easy_pickle.loads(s)
    assert kappa == kappa2

    image = os.path.join(dxtbx_dir, "tests", "omega_scan.cbf")
    assert GoniometerFactory.imgCIF(image)

    kappa = GoniometerFactory.kappa(50.0, -10.0, 30.0, 20.0, "-y", "omega")

    s = easy_pickle.dumps(kappa)
    kappa2 = easy_pickle.loads(s)
    assert kappa == kappa2


def test_multi_axis_goniometer():
    alpha = 50
    omega = -10
    kappa = 30
    phi = 20
    direction = "-y"

    kappa_omega_scan = KappaGoniometer(alpha, omega, kappa, phi, direction, "omega")
    axes = (
        kappa_omega_scan.get_phi_axis(),
        kappa_omega_scan.get_kappa_axis(),
        kappa_omega_scan.get_omega_axis(),
    )
    angles = (
        kappa_omega_scan.get_phi_angle(),
        kappa_omega_scan.get_kappa_angle(),
        kappa_omega_scan.get_omega_angle(),
    )

    # First test a kappa goniometer with omega as the scan axis
    axes = flex.vec3_double(axes)
    angles = flex.double(angles)
    names = flex.std_string(("phi", "kappa", "omega"))
    scan_axis = 2

    multi_axis_omega_scan = GoniometerFactory.multi_axis(axes, angles, names, scan_axis)
    assert approx_equal(
        multi_axis_omega_scan.get_fixed_rotation(),
        kappa_omega_scan.get_fixed_rotation(),
    )
    assert approx_equal(
        multi_axis_omega_scan.get_rotation_axis(), kappa_omega_scan.get_rotation_axis()
    )

    recycle_omega = MultiAxisGoniometer.from_dict(multi_axis_omega_scan.to_dict())
    assert approx_equal(recycle_omega.get_axes(), multi_axis_omega_scan.get_axes())
    assert approx_equal(recycle_omega.get_angles(), multi_axis_omega_scan.get_angles())
    assert recycle_omega.get_scan_axis() == multi_axis_omega_scan.get_scan_axis()

    # Now test a kappa goniometer with phi as the scan axis
    kappa_phi_scan = KappaGoniometer(alpha, omega, kappa, phi, direction, "phi")

    scan_axis = 0
    multi_axis_phi_scan = GoniometerFactory.multi_axis(axes, angles, names, scan_axis)
    assert approx_equal(
        multi_axis_phi_scan.get_fixed_rotation(), kappa_phi_scan.get_fixed_rotation()
    )
    assert approx_equal(
        matrix.sqr(multi_axis_phi_scan.get_setting_rotation())
        * multi_axis_phi_scan.get_rotation_axis_datum(),
        kappa_phi_scan.get_rotation_axis(),
    )
    assert approx_equal(
        multi_axis_phi_scan.get_rotation_axis(), kappa_phi_scan.get_rotation_axis()
    )

    recycle_phi = MultiAxisGoniometer.from_dict(multi_axis_phi_scan.to_dict())
    assert approx_equal(recycle_phi.get_axes(), multi_axis_phi_scan.get_axes())
    assert approx_equal(recycle_phi.get_angles(), multi_axis_phi_scan.get_angles())
    assert recycle_phi.get_scan_axis() == multi_axis_phi_scan.get_scan_axis()

    s = easy_pickle.dumps(multi_axis_phi_scan)
    recycle = easy_pickle.loads(s)
    assert recycle == multi_axis_phi_scan

    assert approx_equal(recycle.get_axes(), multi_axis_phi_scan.get_axes())
    assert approx_equal(recycle.get_angles(), multi_axis_phi_scan.get_angles())
    assert recycle.get_scan_axis() == multi_axis_phi_scan.get_scan_axis()
    recycle.set_angles((0, 90, 180))
    assert approx_equal(recycle.get_angles(), (0, 90, 180))
    new_axes = flex.vec3_double(
        (
            (0.99996, -0.00647, -0.00659),
            (0.91314, 0.27949, -0.29674),
            (1.00000, -0.00013, -0.00064),
        )
    )
    recycle.set_axes(new_axes)
    assert approx_equal(recycle.get_axes(), new_axes)

    # Check exception is raised if scan axis is out range
    try:
        GoniometerFactory.multi_axis(axes, angles, names, 3)
    except RuntimeError:
        pass
    else:
        raise Exception_expected

    # Single axis is just a special case of a multi axis goniometer
    single_axis = GoniometerFactory.multi_axis(
        flex.vec3_double(((1, 0, 0),)), flex.double((0,)), flex.std_string(("PHI",)), 0
    )
    assert single_axis.get_fixed_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)
    assert single_axis.get_setting_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)
    assert single_axis.get_rotation_axis() == (1, 0, 0)


def test_single_axis_goniometer_from_phil():
    params = goniometer_phil_scope.fetch(
        parse(
            """
    goniometer {
      axis = (1, 0, 0)
    }
  """
        )
    ).extract()

    g1 = GoniometerFactory.from_phil(params)

    assert g1.get_rotation_axis() == (1, 0, 0)


def test_single_axis_goniometer_using_axes_from_phil():
    params = goniometer_phil_scope.fetch(
        parse(
            """
    goniometer {
      axes = (1, 0, 0)
    }
  """
        )
    ).extract()

    g1 = GoniometerFactory.from_phil(params)

    assert g1.get_rotation_axis() == (1, 0, 0)


def test_single_axis_goniometer_using_axis_and_axes_from_phil_raises_error():
    """Supplying both the 'axis' and 'axes' parameters is ambiguous, so it should
    be mutually exclusive. This test ensures that is the case."""

    params = goniometer_phil_scope.fetch(
        parse(
            """
    goniometer {
      axis = (1, 0, 0)
      axes = (1, 0, 0)
    }
    """
        )
    ).extract()

    with pytest.raises(ValueError):
        GoniometerFactory.from_phil(params)


def test_single_axis_goniometer_with_fixed_rotation_and_reference_from_phil():
    params = goniometer_phil_scope.fetch(
        parse(
            """
    goniometer {
      axes = (1, 0, 0)
    }
  """
        )
    ).extract()

    g1 = GoniometerFactory.from_phil(params)

    params = goniometer_phil_scope.fetch(
        parse(
            """
    goniometer {
      axis = (0, 1, 0)
      fixed_rotation = (0, 1, 0, 1, 0, 0, 0, 0, 1)
    }
  """
        )
    ).extract()

    g2 = GoniometerFactory.from_phil(params, reference=g1)

    assert g2.get_rotation_axis() == (0, 1, 0)
    assert g2.get_fixed_rotation() == (0, 1, 0, 1, 0, 0, 0, 0, 1)


def test_multi_axis_goniometer_from_phil():
    params = goniometer_phil_scope.fetch(
        parse(
            """
    goniometer {
      axes = (1, 0, 0, 0, 1, 0, 0, 0, 1)
      scan_axis = 2
    }
  """
        )
    ).extract()

    g1 = GoniometerFactory.from_phil(params)
    assert tuple(g1.get_axes()) == ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    assert g1.get_scan_axis() == 2

    # Test using a reference
    params = goniometer_phil_scope.fetch(
        parse(
            """
    goniometer {
      axes = (0, 1, 0, 1, 0, 0, 0, 0, 1)
    }
  """
        )
    ).extract()

    g2 = GoniometerFactory.from_phil(params, reference=g1)

    assert tuple(g2.get_axes()) == ((0, 1, 0), (1, 0, 0), (0, 0, 1))


def test_scan_varying():
    axis = (1, 0, 0)
    g = Goniometer(axis)

    assert g.get_num_scan_points() == 0
    assert g.get_setting_rotation_at_scan_points().size() == 0
    with pytest.raises(RuntimeError):
        g.get_setting_rotation_at_scan_point(0)

    # set varying beam
    num_scan_points = 11
    S_static = matrix.sqr(g.get_setting_rotation())
    S_as_scan_points = [S_static]
    axis = matrix.col.random(3, -1.0, 1.0).normalize()
    R = axis.axis_and_angle_as_r3_rotation_matrix(angle=0.01, deg=True)
    for i in range(num_scan_points - 1):
        S_as_scan_points.append(R * S_as_scan_points[-1])
    g.set_setting_rotation_at_scan_points(S_as_scan_points)
    assert g.get_num_scan_points() == 11
    assert g.get_setting_rotation_at_scan_points().size() == 11

    for t in range(num_scan_points):
        S_t = matrix.sqr(g.get_setting_rotation_at_scan_point(t))
        assert S_t == S_as_scan_points[t]

    # also test setting as tuple
    g.set_setting_rotation_at_scan_points(tuple(S_as_scan_points))
    assert g.get_num_scan_points() == 11
    assert g.get_setting_rotation_at_scan_points().size() == 11

    # test resetting
    g.reset_scan_points()
    assert g.get_num_scan_points() == 0
    assert g.get_setting_rotation_at_scan_points().size() == 0


def test_comparison():
    # Setting rotation for small random offset
    offset_ax = matrix.col.random(3, -1.0, 1.0).normalize()
    S = offset_ax.axis_and_angle_as_r3_rotation_matrix(angle=0.01, deg=True)

    # Equal goniometers with scan-points set
    g1 = Goniometer((1, 0, 0))
    g1.set_setting_rotation(S)
    g1.set_setting_rotation_at_scan_points([S] * 5)
    g2 = Goniometer((1, 0, 0))
    g2.set_setting_rotation(S)
    g2.set_setting_rotation_at_scan_points([S] * 5)

    assert g1 == g2
    assert g1.is_similar_to(g2)

    # Different setting matrix
    g3 = Goniometer((1, 0, 0))
    invS = S.inverse()
    g3.set_setting_rotation(invS)
    g3.set_setting_rotation_at_scan_points([invS] * 5)

    assert g1 != g3
    assert not g1.is_similar_to(g3)
