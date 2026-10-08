from __future__ import annotations

import copy

import pytest

from libtbx.phil import parse
from scitbx import matrix

from dxtbx.model import Beam, BeamBase, PolychromaticBeam, XFELBeam
from dxtbx.model.beam import BeamFactory, Probe, beam_phil_scope


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

    params3 = beam_phil_scope.fetch(
        parse(
            """
    beam {
        probe = electron
        sample_to_source_distance = 4000
    }
  """
        )
    ).extract()
    b4 = BeamFactory.from_phil(params3, reference)
    assert b4.get_probe() == Probe.electron
    assert b4.get_sample_to_source_distance() == pytest.approx(4000)


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


def test_polychromatic_beam_from_phil():
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
      sample_to_source_distance = 5000
      wavelength_range = (0.2, 10)
    }
    """
        )
    ).extract()

    beam = BeamFactory.from_phil(params)
    assert isinstance(beam, PolychromaticBeam)

    assert beam.get_sample_to_source_direction() == pytest.approx((0.0, 0.0, 1.0))
    assert beam.get_divergence() == pytest.approx(0.2)
    assert beam.get_sigma_divergence() == pytest.approx(0.3)
    assert beam.get_polarization_normal() == pytest.approx((0.0, -1.0, 0.0))
    assert beam.get_polarization_fraction() == pytest.approx(0.65)
    assert beam.get_transmission() == pytest.approx(0.5)
    assert beam.get_flux() == pytest.approx(0.75)
    assert beam.get_sample_to_source_distance() == pytest.approx(5000)
    assert beam.get_wavelength_range() == pytest.approx((0.2, 10))


def test_polychromatic_beam_from_dict():
    beam = PolychromaticBeam()
    assert beam == BeamFactory.from_dict(beam.to_dict())


def test_make_polychromatic_beam():
    direction = (0.0, 0.0, 1.0)
    divergence = 0.2
    sigma_divergence = 0.3
    polarization_normal = (0.0, -1.0, 0.0)
    polarization_fraction = 0.65
    transmission = 0.5
    flux = 0.75
    probe = Probe.neutron
    sample_to_source_distance = 8500
    wavelength_range = (0.2, 10)

    beam = BeamFactory.make_polychromatic_beam(
        direction=direction,
        divergence=divergence,
        sigma_divergence=sigma_divergence,
        polarization_normal=polarization_normal,
        polarization_fraction=polarization_fraction,
        transmission=transmission,
        flux=flux,
        probe=probe,
        sample_to_source_distance=sample_to_source_distance,
        wavelength_range=wavelength_range,
    )

    assert beam.get_sample_to_source_direction() == pytest.approx((0.0, 0.0, 1.0))
    assert beam.get_divergence() == pytest.approx(0.2)
    assert beam.get_sigma_divergence() == pytest.approx(0.3)
    assert beam.get_polarization_normal() == pytest.approx((0.0, -1.0, 0.0))
    assert beam.get_polarization_fraction() == pytest.approx(0.65)
    assert beam.get_transmission() == pytest.approx(0.5)
    assert beam.get_flux() == pytest.approx(0.75)
    assert beam.get_probe() == Probe.neutron
    assert beam.get_sample_to_source_distance() == pytest.approx(8500.0)
    assert beam.get_wavelength_range() == pytest.approx((0.2, 10))


def test_polychromatic_beam_wavelength_guards():
    beam = PolychromaticBeam()
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


def test_polychromatic_beam_comparison():
    beam1 = PolychromaticBeam()
    beam1.set_wavelength_range((0.0, 1.0))
    beam1.set_sample_to_source_distance(10)
    beam1.set_direction((0.0, 0.0, 1.0))
    beam2 = PolychromaticBeam()
    beam2.set_wavelength_range((0.0, 1.0))
    beam2.set_sample_to_source_distance(10)
    beam2.set_direction((0.0, 0.0, 1.0))

    assert beam1 == beam2
    beam1.set_wavelength_range((0.0, 2.0))
    assert beam1 != beam2
    beam1.set_wavelength_range((0.0, 1.0))
    assert beam1 == beam2

    beam1.set_sample_to_source_distance(1.0)
    assert beam1 != beam2
    beam1.set_sample_to_source_distance(10)
    assert beam1 == beam2

    beam1.set_probe(Probe.neutron)
    assert beam1 != beam2
    beam1.set_probe(Probe.xray)
    assert beam1 == beam2

    beam1.set_direction((1.0, 0.0, 0.0))
    assert beam1 != beam2
    beam1.set_direction((0.0, 0.0, 1.0))
    assert beam1 == beam2

    beam3 = Beam()
    assert beam1 != beam3


def test_polychromatic_beam_str():
    beam = PolychromaticBeam()
    assert (
        beam.__str__()
        == "Beam:\n    probe: x-ray\n    sample to source direction : {0,0,1}\n    divergence: 0\n    sigma divergence: 0\n    polarization normal: {0,1,0}\n    polarization fraction: 0.5\n    flux: 0\n    transmission: 1\n    sample to source distance : 0\n    wavelength range : {0,0}\n"
    )


def test_copy_beam():
    beam = PolychromaticBeam()
    assert beam == copy.deepcopy(beam)
    beam = Beam()
    assert beam == copy.deepcopy(beam)


def test_beam_hierarchy():
    """The three spectral beam types are direct siblings under BeamBase.

    Beam stays the (only) monochromatic type; PolychromaticBeam and XFELBeam no
    longer inherit from Beam, so generic "any beam" code must bind BeamBase.
    """
    assert issubclass(Beam, BeamBase)
    assert issubclass(PolychromaticBeam, BeamBase)
    assert issubclass(XFELBeam, BeamBase)

    # PolychromaticBeam / XFELBeam are no longer Beam subclasses.
    assert not issubclass(PolychromaticBeam, Beam)
    assert not issubclass(XFELBeam, Beam)

    assert isinstance(PolychromaticBeam(), BeamBase)
    assert isinstance(XFELBeam(), BeamBase)
    assert not isinstance(PolychromaticBeam(), Beam)
    assert not isinstance(XFELBeam(), Beam)

    # The per-shot resolve boundary discriminates by the get_monochromatic_beam
    # capability, so it must live on XFELBeam only -- never on plain Beam.
    assert not hasattr(Beam(), "get_monochromatic_beam")
    assert hasattr(XFELBeam(), "get_monochromatic_beam")


def test_beambase_is_abstract():
    """BeamBase is abstract (no_init) and cannot be instantiated from Python."""
    with pytest.raises((RuntimeError, TypeError)):
        BeamBase()


def test_xfel_beam_spectral_guards():
    """XFELBeam has no fixed wavelength/s0 but reports zero scan points so that
    generic scan-varying probes succeed without XFEL-awareness."""
    beam = XFELBeam()

    # Deliberately 0 (not a throw) -- generic copy/compare/serialize code.
    assert beam.get_num_scan_points() == 0

    with pytest.raises(RuntimeError):
        _ = beam.get_wavelength()
    with pytest.raises(RuntimeError):
        _ = beam.get_s0()
    with pytest.raises(RuntimeError):
        beam.set_wavelength(1.0)
    with pytest.raises(RuntimeError):
        beam.set_s0((0.0, 0.0, 0.1))


def test_xfel_get_monochromatic_beam():
    """get_monochromatic_beam(λ) resolves an XFELBeam to a real monochromatic
    Beam carrying the XFEL geometry/polarization/probe at the supplied λ."""
    beam = XFELBeam(
        direction=(0.0, 0.0, 1.0),
        divergence=0.0,
        sigma_divergence=0.0,
        polarization_normal=(0.0, 1.0, 0.0),
        polarization_fraction=0.6,
        flux=1.5e12,
        transmission=0.8,
        probe=Probe.electron,
        sample_to_source_distance=12.0,
    )

    mono = beam.get_monochromatic_beam(1.3)

    assert isinstance(mono, Beam)
    assert mono.get_wavelength() == pytest.approx(1.3)
    assert mono.get_sample_to_source_direction() == pytest.approx((0.0, 0.0, 1.0))
    assert mono.get_polarization_normal() == pytest.approx((0.0, 1.0, 0.0))
    assert mono.get_polarization_fraction() == pytest.approx(0.6)
    assert mono.get_flux() == pytest.approx(1.5e12)
    assert mono.get_transmission() == pytest.approx(0.8)
    assert mono.get_probe() == Probe.electron
    assert mono.get_sample_to_source_distance() == pytest.approx(12.0)
