from __future__ import annotations

import pytest

from libtbx.phil import parse
from scitbx import matrix

from dxtbx.model import Beam, TOFBeam
from dxtbx.model.beam import (
    BeamBaseFactory,
    BeamFactory,
    BeamType,
    TOFBeamFactory,
    beam_phil_scope,
)


@pytest.fixture
def direction():
    return matrix.col((0.013142, 0.002200, 1.450476))


@pytest.fixture
def wavelength():
    return 0.689400


@pytest.fixture
def sample_to_moderator_distance():
    return 10


@pytest.fixture
def eps():
    return 1e-7


def test_setting_direction_and_wavelength(direction, wavelength, eps):
    unit_direction = direction.normalize()

    # Create the beam
    b = Beam(direction, wavelength)

    # Check direction is a unit vector
    assert matrix.col(b.get_sample_to_source_direction()).length() == pytest.approx(1)
    assert abs(matrix.col(b.get_sample_to_source_direction()) - unit_direction) <= eps

    # Check wavelength is correct
    assert b.get_wavelength() == pytest.approx(wavelength)

    # Check s0 is in direction and has length 1/wavelength
    assert matrix.col(b.get_s0()).length() == pytest.approx(1.0 / wavelength)
    assert abs(-matrix.col(b.get_s0()).normalize() - unit_direction) <= eps


def test_setting_s0(direction, wavelength, eps):
    unit_direction = direction.normalize()
    s0 = -unit_direction * 1.0 / wavelength

    # Create the beam
    b = Beam(s0)

    # Check direction is a unit vector
    assert matrix.col(b.get_sample_to_source_direction()).length() == pytest.approx(1)
    assert abs(matrix.col(b.get_sample_to_source_direction()) - unit_direction) <= eps

    # Check wavelength is correct
    assert b.get_wavelength() == pytest.approx(wavelength)

    # Check s0 is in direction and has length 1/wavelength
    assert matrix.col(b.get_s0()).length() == pytest.approx(1.0 / wavelength)
    assert abs(-matrix.col(b.get_s0()).normalize() - unit_direction) <= eps
    assert abs(matrix.col(b.get_s0()) - s0) <= eps


def test_from_phil(direction, wavelength):

    reference = Beam(direction, wavelength)

    params1 = beam_phil_scope.fetch(
        parse(
            """
    beam {
      wavelength = 1.0
      direction = (0, 0, 1)
    }
  """
        )
    ).extract()

    params2 = beam_phil_scope.fetch(
        parse(
            """
    beam {
      wavelength = 1.0
    }
  """
        )
    ).extract()

    # Create the beam
    assert BeamFactory.from_phil(params1)
    assert BeamFactory.from_phil(params2, reference)
    with pytest.raises(RuntimeError):
        BeamFactory.from_phil(params2)

    params3 = beam_phil_scope.fetch(
        parse(
            """
    beam {
      polarization_normal = 1,0,0
      polarization_fraction = 0.5
    }
  """
        )
    ).extract()
    b3 = BeamFactory.from_phil(params3, reference)
    assert b3.get_polarization_fraction() == 0.5
    assert b3.get_polarization_normal() == (1.0, 0.0, 0.0)


def test_scan_varying(direction, wavelength):
    unit_direction = direction.normalize()
    s0 = -unit_direction * 1.0 / wavelength

    # Create the beam
    b = Beam(s0)

    assert b.get_num_scan_points() == 0
    assert b.get_s0_at_scan_points().size() == 0
    with pytest.raises(RuntimeError):
        b.get_s0_at_scan_point(0)  # should raise RuntimeError

    # set varying beam
    num_scan_points = 11
    s0_static = matrix.col(b.get_s0())
    s0_as_scan_points = [s0_static]
    axis = matrix.col.random(3, -1.0, 1.0).normalize()
    for i in range(num_scan_points - 1):
        s0_as_scan_points.append(
            s0_as_scan_points[-1].rotate_around_origin(axis, angle=0.01, deg=True)
        )
    b.set_s0_at_scan_points(s0_as_scan_points)
    assert b.get_num_scan_points() == 11
    assert b.get_s0_at_scan_points().size() == 11

    for t in range(num_scan_points):
        s0_t = matrix.col(b.get_s0_at_scan_point(t))
        assert s0_t == s0_as_scan_points[t]

    # also test setting as tuple
    b.set_s0_at_scan_points(tuple(s0_as_scan_points))
    assert b.get_num_scan_points() == 11
    assert b.get_s0_at_scan_points().size() == 11

    # test resetting
    b.reset_scan_points()
    assert b.get_num_scan_points() == 0
    assert b.get_s0_at_scan_points().size() == 0


def test_beam_object_comparison(direction, wavelength):
    unit_direction = direction.normalize()
    s0 = -unit_direction * 1.0 / wavelength

    # Equal beams with scan-points set
    b1 = Beam(s0)
    b1.set_s0_at_scan_points([s0] * 5)
    b2 = Beam(s0)
    b2.set_s0_at_scan_points([s0] * 5)

    assert b1 == b2
    assert b1.is_similar_to(b2)

    # Different direction
    b3 = Beam(-s0)
    b3.set_s0_at_scan_points([-s0] * 5)
    assert b1 != b3
    assert not b1.is_similar_to(b3)

    # Different wavelength
    b4 = Beam(s0 * 1.5)
    b4.set_s0_at_scan_points([s0 * 1.5] * 5)
    assert b1 != b4
    assert not b1.is_similar_to(b4)


def test_beam_self_serialization():
    beam = Beam()
    assert beam == BeamFactory.from_dict(beam.to_dict())


def test_beam_factory_make_beam_consistency(
    direction, wavelength, sample_to_moderator_distance
):

    # BeamFactory tests
    beam1 = BeamFactory.make_beam(sample_to_source=direction, wavelength=wavelength)
    beam2 = BeamBaseFactory.make_beam(sample_to_source=direction, wavelength=wavelength)
    beam3 = BeamBaseFactory.make_beam(
        beam_type=BeamType.Monochromatic,
        sample_to_source=direction,
        wavelength=wavelength,
    )

    assert beam1 == beam2 == beam3

    # TOFBeamFactory tests

    beam1 = TOFBeamFactory.make_beam(
        sample_to_source_direction=direction,
        sample_to_moderator_distance=sample_to_moderator_distance,
    )
    beam2 = BeamBaseFactory.make_beam(
        beam_type=BeamType.TOF,
        sample_to_source_direction=direction,
        sample_to_moderator_distance=sample_to_moderator_distance,
    )

    assert beam1 == beam2


def test_beam_factory_from_phil_consistency(
    direction, wavelength, sample_to_moderator_distance
):

    reference = Beam(direction, wavelength)

    # Test BeamBaseFactory defaults to BeamFactory when beam type not present

    params = beam_phil_scope.fetch(
        parse(
            """
    beam {
      wavelength = 1.0
      direction = (0, 0, 1)
    }
    """
        )
    ).extract()

    beam1 = BeamFactory.from_phil(params)
    beam2 = BeamBaseFactory.from_phil(params)
    assert beam1 == beam2

    beam1 = BeamFactory.from_phil(params, reference)
    beam2 = BeamBaseFactory.from_phil(params, reference)
    assert beam1 == beam2

    # Test BeamBaseFactory is consistent with BeamFactory

    params = beam_phil_scope.fetch(
        parse(
            """
    beam {
      type = monochromatic
      wavelength = 1.0
      direction = (0, 0, 1)
    }
    """
        )
    ).extract()

    beam1 = BeamFactory.from_phil(params)
    beam2 = BeamBaseFactory.from_phil(params)
    assert beam1 == beam2

    beam1 = BeamFactory.from_phil(params, reference)
    beam2 = BeamBaseFactory.from_phil(params, reference)
    assert beam1 == beam2

    # Test BeamBaseFactory consistent with TOFBeamFactory

    reference = TOFBeam(direction, sample_to_moderator_distance)

    params = beam_phil_scope.fetch(
        parse(
            """
    beam {
      type = tof
      sample_to_moderator_distance = 1.0
      direction = (0, 0, 1)
    }
    """
        )
    ).extract()

    beam1 = TOFBeamFactory.from_phil(params)
    beam2 = BeamBaseFactory.from_phil(params)
    assert beam1 == beam2

    beam1 = TOFBeamFactory.from_phil(params, reference)
    beam2 = BeamBaseFactory.from_phil(params, reference)
    assert beam1 == beam2


def test_beam_factory_from_dict_consistency():

    beam_dict = Beam().to_dict()
    beam1 = BeamFactory.from_dict(beam_dict)
    beam2 = BeamBaseFactory.from_dict(beam_dict)
    assert beam1 == beam2

    beam_dict = TOFBeam().to_dict()
    beam_dict["sample_to_moderator_distance"] = 10
    beam1 = TOFBeamFactory.from_dict(beam_dict)
    beam2 = BeamBaseFactory.from_dict(beam_dict)
    assert beam1 == beam2


def test_tof_beam_constructor(direction, sample_to_moderator_distance, eps):

    beam = TOFBeam(direction, sample_to_moderator_distance)

    smd = beam.get_sample_to_moderator_distance()
    assert smd == pytest.approx(sample_to_moderator_distance)
    assert matrix.col(beam.get_sample_to_source_direction()).length() == pytest.approx(
        1
    )
    beam_dir = beam.get_sample_to_source_direction()
    assert abs(direction.normalize() - matrix.col(beam_dir)) < eps


def test_tof_beam_factory_make_beam(direction, sample_to_moderator_distance, eps):

    beam = TOFBeamFactory.make_beam(
        sample_to_source_direction=direction,
        sample_to_moderator_distance=sample_to_moderator_distance,
    )

    smd = beam.get_sample_to_moderator_distance()
    assert smd == pytest.approx(sample_to_moderator_distance)
    assert matrix.col(beam.get_sample_to_source_direction()).length() == pytest.approx(
        1
    )
    beam_dir = beam.get_sample_to_source_direction()
    assert abs(direction.normalize() - matrix.col(beam_dir)) < eps


def test_tof_beam_factory_from_phil(eps):

    params = beam_phil_scope.fetch(
        parse(
            """
    beam {
      type = tof
      sample_to_moderator_distance = 10.0
      direction = (1, 2, 3)
    }
    """
        )
    ).extract()

    beam = TOFBeamFactory.from_phil(params)

    expected_dir = matrix.col(params.beam.direction).normalize()
    beam_dir = matrix.col(beam.get_sample_to_source_direction())
    assert abs(beam_dir - expected_dir) < eps
    smd = beam.get_sample_to_moderator_distance()
    assert smd == pytest.approx(params.beam.sample_to_moderator_distance)


def test_tof_beam_factory_from_dict(eps):

    beam_dict = {
        "direction": matrix.col((1, 2, 3)),
        "sample_to_moderator_distance": 10.0,
    }

    beam = TOFBeamFactory.from_dict(dict=beam_dict)

    smd = beam.get_sample_to_moderator_distance()
    assert smd == pytest.approx(beam_dict["sample_to_moderator_distance"])
    assert matrix.col(beam.get_sample_to_source_direction()).length() == pytest.approx(
        1
    )
    beam_dir = beam.get_sample_to_source_direction()
    assert abs(beam_dict["direction"].normalize() - matrix.col(beam_dir)) < eps
