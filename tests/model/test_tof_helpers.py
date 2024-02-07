from __future__ import annotations

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
