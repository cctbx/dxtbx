from __future__ import absolute_import, division, print_function

from builtins import range

import pytest

from libtbx.phil import parse
from scitbx import matrix

from dxtbx.model import Beam, SpectrumBeam
from dxtbx.model.beam import BeamFactory, beam_phil_scope
from dials.array_family import flex


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


def test_spectrum_beam():
    spectrum_energies = flex.double(range(9450, 9550))
    spectrum_weights = flex.double(range(len(spectrum_energies)))
    b1 = SpectrumBeam()
    b2 = SpectrumBeam()
    b1.set_spectrum(spectrum_energies, spectrum_weights)
    b2.set_spectrum(spectrum_energies, spectrum_weights)
    assert b1.get_weighted_wavelength() == pytest.approx(1.3028567060142213)
    assert b1.is_similar_to(b2)
    b2.set_spectrum(spectrum_energies + 50, spectrum_weights)
    assert not b1.is_similar_to(b2)

    b3 = Beam()
    b1.set_wavelength(1.2)
    b3.set_wavelength(1.2)
    assert b1.is_similar_to(b3)
    assert b3.is_similar_to(b1)
