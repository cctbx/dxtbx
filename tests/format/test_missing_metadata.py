"""Tests for Format classes that declare metadata they cannot read from the file."""

from __future__ import annotations

import pytest

from dxtbx.format.Format import Format
from dxtbx.format.FormatStill import FormatStill
from dxtbx.model import MissingMetadataError
from dxtbx.model.beam import BeamFactory
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model.scan import ScanFactory

OVERRIDES = """
geometry {
  beam.wavelength = 0.02508
  goniometer.axis = 1,0,0
  detector.distance = 999
}
"""


class FormatDummy(Format):
    """A format with complete, if fictional, metadata."""

    @staticmethod
    def understand(image_file):
        return False

    def _beam(self):
        return BeamFactory.make_beam(sample_to_source=(0, 0, 1), wavelength=1.0)

    def _detector(self):
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

    def _goniometer(self):
        return GoniometerFactory.known_axis((0, 1, 0))

    def _scan(self):
        return ScanFactory.make_scan(
            (1, 5), 0.1, (0.0, 0.5), dict.fromkeys(range(5), 0)
        )


class FormatDummyED(FormatDummy):
    """The wavelength and rotation axis are not in the headers of this format."""

    missing_metadata = ("beam.wavelength", "goniometer.axis")


class FormatDummyRuntime(FormatDummy):
    """Whether the wavelength is readable depends on the individual file."""

    def _start(self):
        if "noheader" in self._image_file:
            self.declare_missing_metadata("beam.wavelength")


class FormatDummyStill(FormatStill):
    missing_metadata = ("beam.wavelength",)

    @staticmethod
    def understand(image_file):
        return False

    def _beam(self):
        return BeamFactory.make_beam(sample_to_source=(0, 0, 1), wavelength=1.0)

    def _detector(self):
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


class FormatDummyBypassed(FormatDummyED):
    """Rebuilds the beam on every call, discarding anything written to
    self._beam_instance."""

    def get_beam(self, index=None):
        return self._beam()


def test_declaring_nothing_is_unaffected():
    """A format that declares nothing must behave identically, overrides or not."""
    plain = FormatDummy("dummy.img")
    with_overrides = FormatDummy("dummy.img", geometry_overrides=OVERRIDES)
    for fmt in (plain, with_overrides):
        assert fmt.get_beam().get_wavelength() == pytest.approx(1.0)
        assert fmt.get_goniometer().get_rotation_axis() == pytest.approx((0, 1, 0))
        assert fmt.get_detector()[0].get_distance() == pytest.approx(100.0)


def test_missing_metadata_raises_without_overrides():
    with pytest.raises(MissingMetadataError) as excinfo:
        FormatDummyED("dummy.img")
    error = excinfo.value
    assert error.unsatisfied == ("beam.wavelength", "goniometer.axis")
    assert error.format_name == "FormatDummyED"
    assert "dummy.img" in str(error)
    assert "geometry.beam.wavelength=<value>" in str(error)


def test_missing_metadata_raises_for_partial_overrides():
    with pytest.raises(MissingMetadataError) as excinfo:
        FormatDummyED(
            "dummy.img", geometry_overrides="geometry.beam.wavelength=0.02508"
        )
    assert excinfo.value.unsatisfied == ("goniometer.axis",)


def test_overrides_fill_in_declared_metadata():
    fmt = FormatDummyED("dummy.img", geometry_overrides=OVERRIDES)
    assert fmt.get_beam().get_wavelength() == pytest.approx(0.02508)
    assert fmt.get_goniometer().get_rotation_axis() == pytest.approx((1, 0, 0))
    # detector.distance was supplied but not declared, so must be ignored
    assert fmt.get_detector()[0].get_distance() == pytest.approx(100.0)
    # ... as must the undeclared scan
    assert fmt.get_scan().get_oscillation() == pytest.approx((0.0, 0.5))


def test_goniometer_axes_satisfies_declared_axis():
    fmt = FormatDummyED(
        "dummy.img",
        geometry_overrides="""
geometry {
  beam.wavelength = 0.02508
  goniometer.axes = 1,0,0
}
""",
    )
    assert fmt.get_goniometer().get_rotation_axis() == pytest.approx((1, 0, 0))


def test_declare_missing_metadata_at_start_time():
    """A declaration made in _start() only applies to the files that need it."""
    assert FormatDummyRuntime("dummy.img").get_beam().get_wavelength() == pytest.approx(
        1.0
    )
    with pytest.raises(MissingMetadataError):
        FormatDummyRuntime("noheader.img")
    fmt = FormatDummyRuntime("noheader.img", geometry_overrides=OVERRIDES)
    assert fmt.get_beam().get_wavelength() == pytest.approx(0.02508)


def test_declare_missing_metadata_validates():
    fmt = FormatDummy("dummy.img")
    for path in ("detector.panel.pixel_size", "goniometer.invert_rotation_axis"):
        with pytest.raises(ValueError):
            fmt.declare_missing_metadata(path)


def test_still_does_not_swallow_the_error():
    """FormatStill.setup() squashes exceptions, but must not squash this one."""
    with pytest.raises(MissingMetadataError):
        FormatDummyStill("dummy.img")
    fmt = FormatDummyStill("dummy.img", geometry_overrides=OVERRIDES)
    assert fmt.get_beam().get_wavelength() == pytest.approx(0.02508)


def test_bypassed_getter_is_a_loud_error():
    with pytest.raises(TypeError, match="get_beam"):
        FormatDummyBypassed("dummy.img", geometry_overrides=OVERRIDES)


def test_ignore_missing_metadata():
    fmt = FormatDummyED("dummy.img", ignore_missing_metadata=True)
    assert fmt.get_beam().get_wavelength() == pytest.approx(1.0)


class FormatDummyOneImage(FormatDummy):
    """One image per file, as for a multi-file sequence."""

    missing_metadata = ("scan.oscillation",)

    def _scan(self):
        index = int(self._image_file.split("_")[-1].split(".")[0])
        return ScanFactory.make_scan((index, index), 0.1, (0.0, 0.0), {index: 0})


def test_one_image_scans_can_be_summed():
    """Every image of a multi-file sequence gets the same override, so the
    oscillation start must be offset by the image index or the sum fails."""
    overrides = "geometry.scan.oscillation = 10.0,0.5"
    scans = [
        FormatDummyOneImage(f"dummy_{i}.img", geometry_overrides=overrides).get_scan()
        for i in (1, 2, 3)
    ]
    assert scans[0].get_oscillation() == pytest.approx((10.0, 0.5))
    assert scans[1].get_oscillation() == pytest.approx((10.5, 0.5))

    total = scans[0]
    for scan in scans[1:]:
        total += scan
    assert total.get_image_range() == (1, 3)
    assert total.get_oscillation() == pytest.approx((10.0, 0.5))
