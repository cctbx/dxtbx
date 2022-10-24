from __future__ import annotations

import numpy as np
import pytest
from scipy.ndimage import rotate
from scipy.stats import norm

from dxtbx.util import ersatz_uuid4
from dxtbx.util import format_float_with_standard_uncertainty as ffwsu
from dxtbx.util.rotate_and_average import rotate_and_average


def test_format_float_with_standard_uncertainty():
    """Test format_float_with_standard_uncertainty."""

    assert ffwsu(21.234567, 0.0013) == "21.2346(13)"
    assert ffwsu(21.234567, 0.0023) == "21.235(2)"
    assert ffwsu(12345, 45) == "12340(40)"
    assert ffwsu(12.3, 1.2) == "12.3(12)"
    assert ffwsu(-0.2451, 0.8135) == "-0.2(8)"
    assert ffwsu(1.234, 0.196) == "1.2(2)"
    assert ffwsu(1.234, 0.193) == "1.23(19)"
    assert ffwsu(90, 0) == "90"
    assert ffwsu(90.0, 0) == "90.0"
    assert ffwsu(90 + 1e-14, 0) == "90.0"


def test_ersatz_uuid4():
    """Test ersatz UUID4 behaviour:
    - (pseudo)randomness and uniqueness
    - shape e.g. 00000000-0000-0000-0000-000000000000
    - content e.g. hex digits"""

    uuids = [ersatz_uuid4() for j in range(1000)]

    assert len(uuids) == len(set(uuids))

    for uuid in uuids:
        bits = uuid.split("-")
        assert tuple(len(b) for b in bits) == (8, 4, 4, 4, 12)
        for b in bits:
            assert set(b) <= set("0123456789abcdef")


def test_rotate_and_average():
    # create a pseudo-spectrum by multiplying two gaussians, to make
    # final image of 20,200, then rotate by 10 degrees
    x = np.linspace(-50, 49, 200)
    fx = norm.pdf(x, 10, 30)  # mean 10, stddev 30
    y = np.linspace(-10, 9, 20)
    fy = norm.pdf(y, 2, 5)  # mean 2, stddev 5

    xpart = np.vstack([fx] * 20)
    ypart = np.vstack([fy] * 200).transpose()

    img = xpart * ypart  # 2d gaussian
    img = rotate(img, angle=10, reshape=True)

    # unrotate and average along the main axis
    x, y = rotate_and_average(img, -10, deg=True)

    # apply ROI
    roi = slice(15, -15)
    x = x[roi]
    y = y[roi]

    unrotated_mean = np.average(
        range(img.shape[1])[roi], weights=np.mean(img, axis=0)[roi]
    )
    rotated_mean = np.average(x, weights=y)

    assert 110.305114 == pytest.approx(unrotated_mean)
    assert 104.311567 == pytest.approx(rotated_mean)
