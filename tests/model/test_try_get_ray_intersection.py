from __future__ import annotations

import pytest

from scitbx import matrix

from dxtbx.model import Detector, Panel

# ---------------------------------------------------------------------------
# Shared fixture: a tilted panel from a real GXPARM.XDS geometry.
# The same geometry is used in test_ray_intersection.py for the throwing tests.
# ---------------------------------------------------------------------------

_FAST_AXIS = (0.000000, -0.939693, -0.342020)
_SLOW_AXIS = (1.000000, 0.000000, 0.000000)
_NORMAL = (0.000000, -0.342020, 0.939693)
_PIXEL_SIZE = (0.172000, 0.172000)
_PIXEL_ORIGIN = (244.836136, 320.338531)
_IMAGE_SIZE = (487, 619)
_DISTANCE = 122.124901
_WAVELENGTH = 0.689400


def _make_tilted_panel():
    """Return a Panel with real crystallographic geometry."""
    origin = (
        (matrix.col(_FAST_AXIS).normalize() * _PIXEL_SIZE[0]) * (0 - _PIXEL_ORIGIN[0])
        + (matrix.col(_SLOW_AXIS).normalize() * _PIXEL_SIZE[1]) * (0 - _PIXEL_ORIGIN[1])
        + (_DISTANCE * matrix.col(_NORMAL).normalize())
    )
    return Panel(
        "",
        "",
        _FAST_AXIS,
        _SLOW_AXIS,
        origin,
        _PIXEL_SIZE,
        _IMAGE_SIZE,
        (0, 0),
        0.0,
        "",
    )


def _make_simple_panel():
    """Return an axis-aligned panel at z=200 mm, 512×512 pixels."""
    return Panel(
        "",
        "Panel",
        (10, 0, 0),
        (0, 10, 0),
        (0, 0, 200),
        (0.172, 0.172),
        (512, 512),
        (0, 1000),
        0.1,
        "Si",
    )


# ---------------------------------------------------------------------------
# VirtualPanelFrame.try_get_ray_intersection  (inherited by Panel)
# ---------------------------------------------------------------------------


def test_try_get_ray_intersection_hit_returns_same_as_throwing():
    """A ray toward the panel returns the same coord as the throwing variant."""
    panel = _make_tilted_panel()
    # Ray toward the image centre pixel (244.836, 320.339)
    centre_lab = panel.get_pixel_lab_coord(_PIXEL_ORIGIN)
    s1 = (1.0 / _WAVELENGTH) * matrix.col(centre_lab)

    result = panel.try_get_ray_intersection(s1)

    assert result is not None
    expected = panel.get_ray_intersection(s1)
    assert abs(result[0] - expected[0]) < 1e-10
    assert abs(result[1] - expected[1]) < 1e-10


def test_try_get_ray_intersection_hit_values():
    """Hit result matches known mm coordinates for image-centre ray."""
    panel = _make_tilted_panel()
    centre_lab = panel.get_pixel_lab_coord(_PIXEL_ORIGIN)
    s1 = (1.0 / _WAVELENGTH) * matrix.col(centre_lab)

    result = panel.try_get_ray_intersection(s1)

    assert result is not None
    # Known values computed from throwing variant: (42.1118…, 55.0982…)
    assert abs(result[0] - 42.111815392) < 1e-6
    assert abs(result[1] - 55.098227332) < 1e-6


def test_try_get_ray_intersection_miss_returns_none():
    """A ray pointing away from the panel (v[2] <= 0) returns None."""
    panel = _make_tilted_panel()
    origin = matrix.col(panel.get_origin())
    s1_back = -(1.0 / _WAVELENGTH) * origin  # reversed

    result = panel.try_get_ray_intersection(s1_back)

    assert result is None


def test_try_get_ray_intersection_miss_agrees_with_throwing():
    """The same backwards ray that returns None also raises in the throwing variant."""
    panel = _make_tilted_panel()
    origin = matrix.col(panel.get_origin())
    s1_back = -(1.0 / _WAVELENGTH) * origin

    assert panel.try_get_ray_intersection(s1_back) is None
    with pytest.raises(RuntimeError):
        panel.get_ray_intersection(s1_back)


def test_try_get_ray_intersection_no_d_matrix_returns_none():
    """A panel with no D matrix (unset geometry) always returns None."""
    panel = Panel()  # no geometry → D_ is boost::none
    result = panel.try_get_ray_intersection((0, 0, 1))
    assert result is None


# ---------------------------------------------------------------------------
# VirtualPanelFrame.try_get_bidirectional_ray_intersection
# ---------------------------------------------------------------------------


def test_try_get_bidirectional_ray_intersection_hit_forward():
    """Forward ray returns same coord as throwing bidirectional variant."""
    panel = _make_tilted_panel()
    centre_lab = panel.get_pixel_lab_coord(_PIXEL_ORIGIN)
    s1 = (1.0 / _WAVELENGTH) * matrix.col(centre_lab)

    result = panel.try_get_bidirectional_ray_intersection(s1)

    assert result is not None
    expected = panel.get_bidirectional_ray_intersection(s1)
    assert abs(result[0] - expected[0]) < 1e-10
    assert abs(result[1] - expected[1]) < 1e-10


def test_try_get_bidirectional_ray_intersection_backward_ray_hits():
    """Bidirectional allows v[2] < 0, so reversed ray still returns a value."""
    panel = _make_tilted_panel()
    origin = matrix.col(panel.get_origin())
    s1_back = -(1.0 / _WAVELENGTH) * origin

    result = panel.try_get_bidirectional_ray_intersection(s1_back)

    assert result is not None


def test_try_get_bidirectional_ray_intersection_no_d_matrix_returns_none():
    """Panel with no D matrix returns None for bidirectional variant too."""
    panel = Panel()
    result = panel.try_get_bidirectional_ray_intersection((0, 0, 1))
    assert result is None
    with pytest.raises(RuntimeError):
        panel.get_bidirectional_ray_intersection((0, 0, 1))


# ---------------------------------------------------------------------------
# Panel.try_get_ray_intersection_px
# ---------------------------------------------------------------------------


def test_try_get_ray_intersection_px_hit_returns_same_as_throwing():
    """Hit ray: pixel result matches throwing get_ray_intersection_px."""
    panel = _make_tilted_panel()
    centre_lab = panel.get_pixel_lab_coord(_PIXEL_ORIGIN)
    s1 = (1.0 / _WAVELENGTH) * matrix.col(centre_lab)

    result = panel.try_get_ray_intersection_px(s1)

    assert result is not None
    expected = panel.get_ray_intersection_px(s1)
    assert abs(result[0] - expected[0]) < 1e-10
    assert abs(result[1] - expected[1]) < 1e-10


def test_try_get_ray_intersection_px_hit_values():
    """Hit result matches known pixel coordinates for image-centre ray."""
    panel = _make_tilted_panel()
    centre_lab = panel.get_pixel_lab_coord(_PIXEL_ORIGIN)
    s1 = (1.0 / _WAVELENGTH) * matrix.col(centre_lab)

    result = panel.try_get_ray_intersection_px(s1)

    assert result is not None
    # Pixel-space round-trip: should be very close to _PIXEL_ORIGIN
    assert abs(result[0] - _PIXEL_ORIGIN[0]) < 1e-6
    assert abs(result[1] - _PIXEL_ORIGIN[1]) < 1e-6


def test_try_get_ray_intersection_px_miss_returns_none():
    """Backwards ray: pixel variant returns None."""
    panel = _make_tilted_panel()
    origin = matrix.col(panel.get_origin())
    s1_back = -(1.0 / _WAVELENGTH) * origin

    result = panel.try_get_ray_intersection_px(s1_back)

    assert result is None


def test_try_get_ray_intersection_px_miss_agrees_with_throwing():
    """Backwards ray: pixel None result agrees with throwing RuntimeError."""
    panel = _make_tilted_panel()
    origin = matrix.col(panel.get_origin())
    s1_back = -(1.0 / _WAVELENGTH) * origin

    assert panel.try_get_ray_intersection_px(s1_back) is None
    with pytest.raises(RuntimeError):
        panel.get_ray_intersection_px(s1_back)


# ---------------------------------------------------------------------------
# Panel.try_get_bidirectional_ray_intersection_px
# ---------------------------------------------------------------------------


def test_try_get_bidirectional_ray_intersection_px_hit_forward():
    """Forward ray: bidirectional px matches throwing variant."""
    panel = _make_tilted_panel()
    centre_lab = panel.get_pixel_lab_coord(_PIXEL_ORIGIN)
    s1 = (1.0 / _WAVELENGTH) * matrix.col(centre_lab)

    result = panel.try_get_bidirectional_ray_intersection_px(s1)

    assert result is not None
    expected = panel.get_bidirectional_ray_intersection_px(s1)
    assert abs(result[0] - expected[0]) < 1e-10
    assert abs(result[1] - expected[1]) < 1e-10


def test_try_get_bidirectional_ray_intersection_px_backward_ray_hits():
    """Backwards ray: bidirectional px still returns a value (not None)."""
    panel = _make_tilted_panel()
    origin = matrix.col(panel.get_origin())
    s1_back = -(1.0 / _WAVELENGTH) * origin

    result = panel.try_get_bidirectional_ray_intersection_px(s1_back)

    assert result is not None


def test_try_get_bidirectional_ray_intersection_px_no_d_matrix_returns_none():
    """Panel with no D matrix returns None for bidirectional px variant."""
    panel = Panel()
    result = panel.try_get_bidirectional_ray_intersection_px((0, 0, 1))
    assert result is None
    with pytest.raises(RuntimeError):
        panel.get_bidirectional_ray_intersection_px((0, 0, 1))


# ---------------------------------------------------------------------------
# Detector.try_get_ray_intersection  (single-ray overload)
# ---------------------------------------------------------------------------


def test_detector_try_get_ray_intersection_hit_returns_panel_and_coord():
    """Hit ray: result is (panel_index, mm_coord), panel_index is 0."""
    panel = _make_simple_panel()
    detector = Detector(panel)
    # Ray along z-axis hits the panel at z=200, landing at mm (0, 0) origin.
    s1 = (0, 0, 1)

    result = detector.try_get_ray_intersection(s1)

    assert result is not None
    panel_idx, mm_coord = result
    assert panel_idx == 0
    assert abs(mm_coord[0]) < 1e-10
    assert abs(mm_coord[1]) < 1e-10


def test_detector_try_get_ray_intersection_hit_matches_throwing():
    """Hit ray: try result matches the throwing get_ray_intersection."""
    panel = _make_simple_panel()
    detector = Detector(panel)
    s1 = (0, 0, 1)

    result = detector.try_get_ray_intersection(s1)
    expected = detector.get_ray_intersection(s1)

    assert result is not None
    assert result[0] == expected[0]
    assert abs(result[1][0] - expected[1][0]) < 1e-10
    assert abs(result[1][1] - expected[1][1]) < 1e-10


def test_detector_try_get_ray_intersection_miss_returns_none():
    """Ray pointing away from all panels returns None."""
    panel = _make_simple_panel()
    detector = Detector(panel)
    s1 = (0, 0, -1)  # points in -z, panel is at +z

    result = detector.try_get_ray_intersection(s1)

    assert result is None


def test_detector_try_get_ray_intersection_miss_agrees_with_throwing():
    """Missed ray: None result agrees with RuntimeError from throwing variant."""
    panel = _make_simple_panel()
    detector = Detector(panel)
    s1 = (0, 0, -1)

    assert detector.try_get_ray_intersection(s1) is None
    with pytest.raises(RuntimeError):
        detector.get_ray_intersection(s1)
