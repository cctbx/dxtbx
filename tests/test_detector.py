import os

import libtbx.load_env

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

    dxtbx_dir = libtbx.env.dist_path("dxtbx")

    image = os.path.join(dxtbx_dir, "tests", "phi_scan_001.cbf")
    # xparm = os.path.join(dxtbx_dir, "tests", "example-xparm.xds")

    assert DetectorFactory.imgCIF(image, "CCD")
    # x = DetectorFactory.XDS(xparm)
