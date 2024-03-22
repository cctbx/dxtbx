from __future__ import annotations

from typing import List

from scipy.constants import Planck, m_n
from scipy.interpolate import interp1d

import cctbx.array_family.flex as flex


def wavelength_from_tof(
    distance: float | flex.double, tof: float | flex.double
) -> float | flex.double:
    """
    distance (m)
    tof (s)

    return (A)
    """
    return ((Planck * tof) / (m_n * distance)) * 10**10


def tof_from_wavelength(
    distance: float | flex.double, wavelength: float | flex.double
) -> float | flex.double:
    """
    wavelength (A)
    return (s)
    """

    wavelength = wavelength * 10**-10

    return (wavelength * m_n * distance) / Planck


def frame_to_tof_interpolator(frames: List[float], tof: List[float]) -> interp1d:

    """
    ToF can vary nonlinearly with frame number.
    A cubic spline is used to better account for these cases.
    """
    assert min(frames) >= 0
    assert min(tof) >= 0
    assert len(frames) == len(tof)
    assert all(i < j for i, j in zip(frames, frames[1:]))
    assert all(i < j for i, j in zip(tof, tof[1:]))
    return interp1d(frames, tof, kind="cubic")


def tof_to_frame_interpolator(tof: List[float], frames: List[float]) -> interp1d:

    """
    ToF can vary nonlinearly with frame number.
    A cubic spline is used to better account for these cases.
    """
    assert min(frames) >= 0
    assert min(tof) >= 0
    assert len(frames) == len(tof)
    assert all(i < j for i, j in zip(frames, frames[1:]))
    assert all(i < j for i, j in zip(tof, tof[1:]))
    return interp1d(tof, frames, kind="cubic")
