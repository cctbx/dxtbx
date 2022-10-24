from __future__ import annotations

from pathlib import Path

from dxtbx.model.detector import DetectorFactory


def test_detector():
    """A test class for the detector class."""

    assert DetectorFactory.simple(
        "CCD",
        100.0,
        (45.0, 52.0),
        "+x",
        "-y",
        (0.172, 0.172),
        (516, 590),
        (0, 1024),
        [],
    )
    assert DetectorFactory.two_theta(
        "CCD",
        60.0,
        (35.0, 34.0),
        "+x",
        "+y",
        "+x",
        30,
        (0.07, 0.07),
        (1042, 1042),
        (0, 1024),
        [],
    )

    image = Path(__file__).parent / "phi_scan_001.cbf"

    assert DetectorFactory.imgCIF(str(image), "CCD")
    # x = DetectorFactory.XDS(xparm)
