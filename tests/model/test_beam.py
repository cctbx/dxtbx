from __future__ import annotations

import pytest

from libtbx.phil import parse
from scitbx import matrix

from dxtbx.model import Beam, PolyBeam
from dxtbx.model.beam import BeamFactory, beam_phil_scope


def test_setting_direction_and_wavelength():
    direction = matrix.col((0.013142, 0.002200, 1.450476))
    unit_direction = direction.normalize()
    wavelength = 0.689400

    # Create the beam
    b = Beam(direction, wavelength)

    eps = 1e-7

    # Check direction is a unit vector
    assert matrix.col(b.get_sample_to_source_direction()).length() == pytest.approx(1)
    assert abs(matrix.col(b.get_sample_to_source_direction()) - unit_direction) <= eps

    # Check wavelength is correct
    assert b.get_wavelength() == pytest.approx(wavelength)

    # Check s0 is in direction and has length 1/wavelength
    assert matrix.col(b.get_s0()).length() == pytest.approx(1.0 / wavelength)
    assert abs(-matrix.col(b.get_s0()).normalize() - unit_direction) <= eps


def test_setting_s0():
    direction = matrix.col((0.013142, 0.002200, 1.450476))
    unit_direction = direction.normalize()
    wavelength = 0.689400
    s0 = -unit_direction * 1.0 / wavelength

    # Create the beam
    b = Beam(s0)

    eps = 1e-7

    # Check direction is a unit vector
    assert matrix.col(b.get_sample_to_source_direction()).length() == pytest.approx(1)
    assert abs(matrix.col(b.get_sample_to_source_direction()) - unit_direction) <= eps

    # Check wavelength is correct
    assert b.get_wavelength() == pytest.approx(wavelength)

    # Check s0 is in direction and has length 1/wavelength
    assert matrix.col(b.get_s0()).length() == pytest.approx(1.0 / wavelength)
    assert abs(-matrix.col(b.get_s0()).normalize() - unit_direction) <= eps
    assert abs(matrix.col(b.get_s0()) - s0) <= eps


def test_from_phil():
    direction = matrix.col((0.013142, 0.002200, 1.450476))
    wavelength = 0.689400

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


def test_scan_varying():
    direction = matrix.col((0.013142, 0.002200, 1.450476))
    unit_direction = direction.normalize()
    wavelength = 0.689400
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


def test_beam_object_comparison():
    direction = matrix.col((0.013142, 0.002200, 1.450476))
    unit_direction = direction.normalize()
    wavelength = 0.689400
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


def test_polybeam_from_phil():
    params = beam_phil_scope.fetch(
        parse(
            """
    beam {
      type = polychromatic
      direction = (0., 0., 1.)
      divergence = 0.2
      sigma_divergence = 0.3
      polarization_normal = (0., -1., 0.)
      polarization_fraction = .65
      transmission = .5
      flux = .75
    }
    """
        )
    ).extract()

    beam = BeamFactory.from_phil(params)
    assert isinstance(beam, PolyBeam)

    assert beam.get_sample_to_source_direction() == pytest.approx((0.0, 0.0, 1.0))
    assert beam.get_divergence() == pytest.approx(0.2)
    assert beam.get_sigma_divergence() == pytest.approx(0.3)
    assert beam.get_polarization_normal() == pytest.approx((0.0, -1.0, 0.0))
    assert beam.get_polarization_fraction() == pytest.approx(0.65)
    assert beam.get_transmission() == pytest.approx(0.5)
    assert beam.get_flux() == pytest.approx(0.75)


def test_polybeam_from_dict():
    beam = PolyBeam()
    assert beam == BeamFactory.from_dict(beam.to_dict())


def test_make_polybeam():

    direction = (0.0, 0.0, 1.0)
    divergence = 0.2
    sigma_divergence = 0.3
    polarization_normal = (0.0, -1.0, 0.0)
    polarization_fraction = 0.65
    transmission = 0.5
    flux = 0.75

    beam = BeamFactory.make_polybeam(
        direction=direction,
        divergence=divergence,
        sigma_divergence=sigma_divergence,
        polarization_normal=polarization_normal,
        polarization_fraction=polarization_fraction,
        transmission=transmission,
        flux=flux,
    )

    assert beam.get_sample_to_source_direction() == pytest.approx((0.0, 0.0, 1.0))
    assert beam.get_divergence() == pytest.approx(0.2)
    assert beam.get_sigma_divergence() == pytest.approx(0.3)
    assert beam.get_polarization_normal() == pytest.approx((0.0, -1.0, 0.0))
    assert beam.get_polarization_fraction() == pytest.approx(0.65)
    assert beam.get_transmission() == pytest.approx(0.5)
    assert beam.get_flux() == pytest.approx(0.75)


def test_polybeam_wavelength_guards():
    beam = PolyBeam()
    with pytest.raises(RuntimeError):
        _ = beam.get_wavelength()
    with pytest.raises(RuntimeError):
        _ = beam.get_s0()
    with pytest.raises(RuntimeError):
        _ = beam.get_num_scan_points()
    with pytest.raises(RuntimeError):
        _ = beam.get_s0_at_scan_points()
    with pytest.raises(RuntimeError):
        _ = beam.get_s0_at_scan_point(0)
    with pytest.raises(RuntimeError):
        beam.reset_scan_points()
    with pytest.raises(RuntimeError):
        beam.set_wavelength(1.0)
    with pytest.raises(RuntimeError):
        beam.set_s0((0.0, 0.0, 0.1))


def test_polybeam_str():
    beam = PolyBeam()
    assert (
        beam.__str__()
        == "Beam:\n    sample to source direction : {0,0,1}\n    divergence: 0\n    sigma divergence: 0\n    polarization normal: {0,1,0}\n    polarization fraction: 0.999\n    flux: 0\n    transmission: 1\n    sample to source distance (m): 0\n"
    )
