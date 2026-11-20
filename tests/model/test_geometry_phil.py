from __future__ import annotations

import pytest

import libtbx.phil

from dxtbx.model import MissingMetadataError, geometry_phil_scope
from dxtbx.model.beam import BeamFactory, beam_phil_scope
from dxtbx.model.detector import DetectorFactory, detector_phil_scope
from dxtbx.model.geometry import (
    apply_geometry_overrides,
    declarable_metadata_paths,
    geometry_params_for_paths,
    select_geometry_phil,
    validate_metadata_path,
)
from dxtbx.model.goniometer import goniometer_phil_scope
from dxtbx.model.scan import scan_phil_scope


def make_detector():
    return DetectorFactory.simple(
        "PAD",
        100.0,
        (50.0, 50.0),
        "+x",
        "-y",
        (0.1, 0.1),
        (500, 500),
        (0, 200000),
        [],
    )


def test_fetch_empties_multiple_scopes():
    """The load-bearing invariant: only fetch() empties the repeatable scopes.

    geometry_phil_scope.extract() materialises one bogus detector.panel and one
    hierarchy.group with id=None, which makes DetectorFactory.from_phil fail.
    Everything in this module must therefore go through fetch().
    """
    bad = geometry_phil_scope.extract().geometry
    assert len(bad.detector.panel) == 1
    assert len(bad.detector.hierarchy.group) == 1
    assert bad.detector.hierarchy.group[0].id is None

    empty = libtbx.phil.scope(name="", objects=[])
    good = geometry_phil_scope.fetch(source=empty).extract().geometry
    assert good.detector.panel == []
    assert good.detector.hierarchy.group == []

    # ... and the extract that fails is precisely the one the factory chokes on
    with pytest.raises(Exception):
        DetectorFactory.from_phil(bad, make_detector())
    assert DetectorFactory.from_phil(good, make_detector()) is not None


def test_geometry_phil_scope_matches_component_scopes():
    composed = {definition.path for definition in geometry_phil_scope.all_definitions()}
    expected = set()
    for scope in (
        beam_phil_scope,
        detector_phil_scope,
        goniometer_phil_scope,
        scan_phil_scope,
    ):
        expected.update(
            f"geometry.{definition.path}" for definition in scope.all_definitions()
        )
    assert composed == expected


def test_declarable_metadata_paths_excludes_modifiers():
    paths = declarable_metadata_paths()
    assert "beam.wavelength" in paths
    assert "detector.distance" in paths
    assert "goniometer.axis" in paths
    assert "scan.oscillation" in paths
    assert "goniometer.invert_rotation_axis" not in paths
    assert "scan.extrapolate_scan" not in paths
    assert "scan.image_range" not in paths
    assert not any(".panel" in path or ".hierarchy" in path for path in paths)


@pytest.mark.parametrize(
    "path",
    [
        "detector.panel.pixel_size",
        "detector.hierarchy.group.id",
        "goniometer.invert_rotation_axis",
        "scan.image_range",
        "beam.not_a_parameter",
        "geometry.beam.wavelength",  # must be relative to the geometry scope
    ],
)
def test_validate_metadata_path_rejects(path):
    with pytest.raises(ValueError):
        validate_metadata_path(path)


def test_validate_metadata_path_accepts():
    for path in ("beam.wavelength", "goniometer.axes", "detector.distance"):
        validate_metadata_path(path)


def test_select_geometry_phil_masks():
    source = libtbx.phil.parse(
        """
geometry {
  beam.wavelength = 0.02508
  detector.distance = 999
}
"""
    )
    masked, supplied = select_geometry_phil(source, ["beam.wavelength"])
    assert supplied == {"beam.wavelength"}
    paths = {definition.path for definition in masked.all_definitions()}
    assert paths == {"geometry.beam.wavelength"}


def test_geometry_params_for_paths():
    overrides = """
geometry {
  beam.wavelength = 0.02508
  detector.distance = 999
}
"""
    params, supplied = geometry_params_for_paths(overrides, ["beam.wavelength"])
    assert supplied == {"beam.wavelength"}
    assert params.beam.wavelength == pytest.approx(0.02508)
    assert params.detector.distance is None
    assert params.detector.panel == []


def test_geometry_params_for_paths_with_no_overrides():
    params, supplied = geometry_params_for_paths(None, ["beam.wavelength"])
    assert supplied == set()
    assert params.beam.wavelength is None
    assert params.detector.panel == []


def test_goniometer_axis_alias():
    for supplied_name in ("axis", "axes"):
        overrides = f"geometry.goniometer.{supplied_name} = 1,0,0"
        for declared in ("goniometer.axis", "goniometer.axes"):
            _, supplied = geometry_params_for_paths(overrides, [declared])
            assert supplied == {f"goniometer.{supplied_name}"}


def test_overwrite_from_phil_is_noop_with_masked_params():
    """A masked extract must not disturb parts of the detector it says nothing about."""
    from dxtbx.model import ParallaxCorrectedPxMmStrategy

    detector = make_detector()
    detector[0].set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(0.5, 0.32))
    origin = detector[0].get_origin()

    params, _ = geometry_params_for_paths(
        "geometry.beam.wavelength = 0.02508", ["beam.wavelength"]
    )
    DetectorFactory.from_phil(params, detector)

    assert isinstance(detector[0].get_px_mm_strategy(), ParallaxCorrectedPxMmStrategy)
    assert detector[0].get_origin() == pytest.approx(origin)


def test_apply_geometry_overrides_only_touches_declared():
    beam = BeamFactory.make_beam(sample_to_source=(0, 0, 1), wavelength=1.0)
    detector = make_detector()
    overrides = """
geometry {
  beam.wavelength = 0.02508
  detector.distance = 999
}
"""
    beam, detector, _, _ = apply_geometry_overrides(
        {"beam.wavelength"}, overrides, beam=beam, detector=detector
    )
    assert beam.get_wavelength() == pytest.approx(0.02508)
    assert detector[0].get_distance() == pytest.approx(100.0)


def test_apply_geometry_overrides_no_declaration_is_identity():
    beam = BeamFactory.make_beam(sample_to_source=(0, 0, 1), wavelength=1.0)
    detector = make_detector()
    result = apply_geometry_overrides(
        set(), "geometry.beam.wavelength = 0.02508", beam=beam, detector=detector
    )
    assert result == (beam, detector, None, None)
    assert beam.get_wavelength() == pytest.approx(1.0)


def test_apply_geometry_overrides_raises_for_unsatisfied():
    beam = BeamFactory.make_beam(sample_to_source=(0, 0, 1), wavelength=1.0)
    with pytest.raises(MissingMetadataError) as excinfo:
        apply_geometry_overrides(
            {"beam.wavelength", "goniometer.axis"},
            "geometry.beam.wavelength = 0.02508",
            beam=beam,
            format_name="FormatDummy",
            image_file="dummy.img",
        )
    error = excinfo.value
    assert error.unsatisfied == ("goniometer.axis",)
    assert "FormatDummy" in str(error)
    assert "dummy.img" in str(error)
    assert "geometry.goniometer.axis=<value>" in str(error)
    # only the unsatisfied path is reported
    assert "geometry.beam.wavelength" not in str(error)


def test_apply_geometry_overrides_none_beam_needs_direction():
    with pytest.raises(MissingMetadataError) as excinfo:
        apply_geometry_overrides(
            {"beam.wavelength"},
            "geometry.beam.wavelength = 0.02508",
            beam=None,
            format_name="FormatDummy",
        )
    assert "beam.direction" in str(excinfo.value)


def test_apply_geometry_overrides_none_detector_is_rejected():
    with pytest.raises(MissingMetadataError) as excinfo:
        apply_geometry_overrides(
            {"detector.distance"},
            "geometry.detector.distance = 999",
            detector=None,
            format_name="FormatDummy",
        )
    assert "placeholder detector" in str(excinfo.value)
