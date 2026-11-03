from __future__ import annotations

import numpy as np
import pytest

from dxtbx.model import tof_helpers


def test_wavelength_from_tof():
    L0 = 8.3
    L1 = 0.24062553870557318
    tof = 0.0017464535727308615
    wavelength = tof_helpers.wavelength_from_tof(L0 + L1, tof)
    assert wavelength == pytest.approx(0.8089606203802903)

    L1 = 0.23021894342637692
    tof = 0.0016982786984656873
    wavelength = tof_helpers.wavelength_from_tof(L0 + L1, tof)
    assert wavelength == pytest.approx(0.7876056098518008)


def test_tof_from_wavelength():
    L0 = 8.3
    L1 = 0.24062553870557318
    wavelength = 0.8089606203802903
    tof = tof_helpers.tof_from_wavelength(L0 + L1, wavelength)
    assert tof == pytest.approx(0.0017464535727308615)

    L1 = 0.23021894342637692
    wavelength = 0.7876056098518008
    tof = tof_helpers.tof_from_wavelength(L0 + L1, wavelength)
    assert tof == pytest.approx(0.0016982786984656873)


def test_frame_to_tof_interpolator():
    frames = [i + 1 for i in range(9)]
    tof = [
        500.5,
        501.5010070800781,
        502.5040283203125,
        503.5090026855469,
        504.5159912109375,
        505.5249938964844,
        506.5360107421875,
        507.54901123046875,
        508.5639953613281,
    ]
    interpolator = tof_helpers.frame_to_tof_interpolator(frames, tof)

    assert interpolator(1) == pytest.approx(500.5)
    assert interpolator(9) == pytest.approx(508.56399536)
    assert interpolator(5.5) == pytest.approx(505.02024061)
    assert interpolator([1, 9, 5.5]) == pytest.approx(
        (500.5, 508.56399536, 505.02024061)
    )
    with pytest.raises(ValueError):
        _ = interpolator(-1)
    with pytest.raises(ValueError):
        _ = interpolator(0)
    with pytest.raises(ValueError):
        _ = interpolator(10)
    with pytest.raises(ValueError):
        _ = interpolator([-1, 2, 3])

    with pytest.raises(AssertionError):
        frames2 = frames + [10]  # frames and tof different length
        _ = tof_helpers.frame_to_tof_interpolator(frames2, tof)
    with pytest.raises(AssertionError):
        frames2 = frames[:]
        frames2[-1] = 5  # frames not strictly increasing
        _ = tof_helpers.frame_to_tof_interpolator(frames2, tof)
    with pytest.raises(AssertionError):
        frames2 = frames[:]  # frames not all > 0
        frames2[0] = -1
        _ = tof_helpers.frame_to_tof_interpolator(frames2, tof)
    with pytest.raises(AssertionError):
        tof2 = tof[:]
        tof2[-1] = 5  # tof not strictly increasing
        _ = tof_helpers.frame_to_tof_interpolator(frames, tof2)
    with pytest.raises(AssertionError):
        tof2 = tof[:]  # tof not all > 0
        tof2[0] = -1
        _ = tof_helpers.frame_to_tof_interpolator(frames, tof2)


def test_tof_to_frame_interpolator():
    frames = [i + 1 for i in range(9)]
    tof = [
        500.5,
        501.5010070800781,
        502.5040283203125,
        503.5090026855469,
        504.5159912109375,
        505.5249938964844,
        506.5360107421875,
        507.54901123046875,
        508.5639953613281,
    ]
    interpolator = tof_helpers.tof_to_frame_interpolator(tof, frames)

    assert interpolator(500.5) == pytest.approx(1)
    assert interpolator(508.56399536) == pytest.approx(9)
    assert interpolator(505.02024061) == pytest.approx(5.5)
    assert interpolator((500.5, 508.56399536, 505.02024061)) == pytest.approx(
        (1, 9, 5.5)
    )
    with pytest.raises(ValueError):
        _ = interpolator(-1)
    with pytest.raises(ValueError):
        _ = interpolator(500)
    with pytest.raises(ValueError):
        _ = interpolator(510)
    with pytest.raises(ValueError):
        _ = interpolator([-1, 500.5, 503])

    with pytest.raises(AssertionError):
        frames2 = frames + [10]  # frames and tof different length
        _ = tof_helpers.tof_to_frame_interpolator(tof, frames2)
    with pytest.raises(AssertionError):
        frames2 = frames[:]
        frames2[-1] = 5  # frames not strictly increasing
        _ = tof_helpers.tof_to_frame_interpolator(tof, frames2)
    with pytest.raises(AssertionError):
        frames2 = frames[:]  # frames not all > 0
        frames2[0] = -1
        _ = tof_helpers.tof_to_frame_interpolator(tof, frames2)
    with pytest.raises(AssertionError):
        tof2 = tof[:]
        tof2[-1] = 5  # tof not strictly increasing
        _ = tof_helpers.tof_to_frame_interpolator(tof2, frames)
    with pytest.raises(AssertionError):
        tof2 = tof[:]  # tof not all > 0
        tof2[0] = -1
        _ = tof_helpers.tof_to_frame_interpolator(tof2, frames)


DETCAL = """\
# NEW CALIBRATION FILE FORMAT (in NeXus/SNS coordinates):
# Lengths are in centimeters.
6         L1     T0_SHIFT
7  3007.0476        0.000
4 DETNUM  NROWS  NCOLS   WIDTH   HEIGHT   DEPTH   DETD   CenterX   CenterY   CenterZ    BaseX    BaseY    BaseZ      UpX      UpY      UpZ
5      1    256    256 15.8208 15.8208   0.2000  40.75    0.9822  -37.2112   16.5912 -0.99976 -0.02115  0.00561 -0.00185  0.33725  0.94141
5     28    256    256 15.8411 15.8727   0.2000  42.59  -32.8706  -12.5964   23.9787 -0.59741  0.66233 -0.45212  0.27120  0.69743  0.66335
"""


def _make_imageset(panel_axes):
    from dxtbx.format.Format import Format
    from dxtbx.imageset import ImageSetFactory
    from dxtbx.model import BeamFactory, Detector

    detector = Detector()
    root = detector.hierarchy()
    for name, (fast, slow) in panel_axes.items():
        panel = root.add_group().add_panel()
        panel.set_name(name)
        panel.set_image_size((256, 256))
        panel.set_pixel_size((0.618, 0.618))
        fast = np.array(fast) / np.linalg.norm(fast)
        slow = np.array(slow) - np.dot(slow, fast) * fast
        slow /= np.linalg.norm(slow)
        panel.set_local_frame(tuple(fast), tuple(slow), (0.0, 0.0, 400.0))
    beam = BeamFactory.make_polychromatic_beam(
        direction=(0, 0, -1),
        wavelength_range=(2.0, 4.0),
        sample_to_source_distance=30000,
    )
    imageset = ImageSetFactory.make_imageset(
        ["dummy.nxs"], format_class=Format, check_format=False
    )
    imageset.set_detector(detector)
    imageset.set_beam(beam)
    return imageset


def _nominal_axes():
    return {
        "bank1": ((-1, 0, 0), (0, 0.34, 0.94)),
        "bank28": ((-0.6, 0.66, -0.45), (0.27, 0.70, 0.66)),
        "bank99": ((1, 0, 0), (0, 1, 0)),
    }


def test_read_detcal(tmp_path):
    detcal = tmp_path / "test.DetCal"
    detcal.write_text(DETCAL)
    instrument, panels = tof_helpers.read_detcal(str(detcal))
    assert instrument == {"L1": 3007.0476, "T0_SHIFT": 0.0}
    assert sorted(panels) == [1, 28]
    assert panels[28]["WIDTH"] == 15.8411
    assert panels[28]["HEIGHT"] == 15.8727
    assert panels[1]["CENTERZ"] == 16.5912


def test_read_detcal_errors(tmp_path):
    detcal = tmp_path / "test.DetCal"
    detcal.write_text(DETCAL.replace("0.94141\n", "\n"))
    with pytest.raises(ValueError, match="expected 16 values, got 15"):
        tof_helpers.read_detcal(detcal)

    detcal.write_text("\n".join(line for line in DETCAL.splitlines() if line[0] != "4"))
    with pytest.raises(ValueError, match="No panel records"):
        tof_helpers.read_detcal(detcal)


def test_apply_detcal_geometry(tmp_path):
    detcal = tmp_path / "test.DetCal"
    detcal.write_text(DETCAL)
    imageset = _make_imageset(_nominal_axes())
    tof_helpers.apply_detcal_geometry(imageset, str(detcal))

    detector = imageset.get_detector()
    bank28 = detector[1]
    assert bank28.get_pixel_size() == pytest.approx((158.411 / 256, 158.727 / 256))
    # The panel centre lands on the calibrated Center
    w, h = 158.411, 158.727
    centre = (
        np.array(bank28.get_origin())
        + 0.5 * w * np.array(bank28.get_fast_axis())
        + 0.5 * h * np.array(bank28.get_slow_axis())
    )
    assert centre == pytest.approx([-328.706, -125.964, 239.787])
    assert np.dot(bank28.get_fast_axis(), bank28.get_slow_axis()) == pytest.approx(0)

    # Panels without a DetCal entry are left alone
    assert detector[2].get_origin() == pytest.approx((0, 0, 400))
    assert imageset.get_beam().get_sample_to_source_distance() == pytest.approx(
        30070.476
    )


def test_apply_detcal_geometry_errors(tmp_path):
    detcal = tmp_path / "test.DetCal"
    detcal.write_text(DETCAL)

    axes = _nominal_axes()
    axes["bank1"] = ((1, 0, 0), (0, 1, 0))
    imageset = _make_imageset(axes)
    with pytest.raises(ValueError, match="disagree by more than 25 deg"):
        tof_helpers.apply_detcal_geometry(imageset, detcal)
    # Nothing was changed, including the panel that did match
    assert imageset.get_detector()[1].get_origin() == pytest.approx((0, 0, 400))

    imageset = _make_imageset({"bank99": ((1, 0, 0), (0, 1, 0))})
    with pytest.raises(ValueError, match="matched none"):
        tof_helpers.apply_detcal_geometry(imageset, detcal)

    detcal.write_text(DETCAL.replace("40.75", "45.75"))
    with pytest.raises(ValueError, match="DETD"):
        tof_helpers.apply_detcal_geometry(_make_imageset(_nominal_axes()), detcal)

    detcal.write_text(
        DETCAL.replace("3007.0476        0.000", "3007.0476        1.500")
    )
    with pytest.raises(NotImplementedError, match="T0_SHIFT"):
        tof_helpers.apply_detcal_geometry(_make_imageset(_nominal_axes()), detcal)
