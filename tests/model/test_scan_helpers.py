from __future__ import annotations

import random
from math import isclose

from dxtbx.model import (
    get_mod2pi_angles_in_range,
    get_range_of_mod2pi_angles,
    is_angle_in_range,
)


def test_is_angle_in_random_range():
    """Test that for a range of angles and angular ranges, the
    is_angle_in_range function correctly calculates if the angle
    is in the range."""

    # Create a number of ranges between 0 and 360 to see if they are within range
    num_range = 100
    total_range = list(range(-720, 720))

    for n in range(num_range):

        # Possible range between (-360, 360)
        angular_range = [
            int(random.random() * 360) - int(random.random() * 360),
            int(random.random() * 360) - int(random.random() * 360),
        ]

        mod_angular_range = [angular_range[0] % 360, angular_range[1] % 360]
        if angular_range[0] < 0 and isclose(angular_range[1], 0):
            mod_angular_range[1] = 360
        if angular_range[1] < 0 and isclose(angular_range[0], 0):
            mod_angular_range[0] = 360

        range_start = min(mod_angular_range[0], mod_angular_range[1])
        range_end = max(mod_angular_range[0], mod_angular_range[1])

        expected_range = list(range(range_start, range_end + 1))

        # Special case for 360, where 0 is also valid
        if 360 in expected_range:
            expected_range.append(0)

        for i in total_range:
            if i % 360 in expected_range:
                assert is_angle_in_range(mod_angular_range, i % 360, True)
            else:
                assert is_angle_in_range(mod_angular_range, i % 360, True) is False


def test_is_angle_sequence_in_range():
    # Create a range over 360 and make sure all angles are valid
    angular_range = (-10, 370)
    for angle in range(0, 360):
        assert is_angle_in_range(angular_range, angle, True) is True


def test_get_range_of_mod2pi_angles():
    """Get the range of equivalent within a given angular range."""
    # In a 360 deg range, have only 1 angle
    a0, a1 = get_range_of_mod2pi_angles((0, 360), 180, deg=True)
    assert a0 == 180 and a1 == 180

    # If no angles within range, a0 > a1
    a0, a1 = get_range_of_mod2pi_angles((181, 360), 180, deg=True)
    assert a0 > a1

    # With 720 deg range, have 2 angles
    a0, a1 = get_range_of_mod2pi_angles((0, 720), 180, deg=True)
    assert a0 == 180 and a1 == 540

    # With 1080 deg range, have 3 angles
    a0, a1 = get_range_of_mod2pi_angles((-360, 720), 180, deg=True)
    assert a0 == -180 and a1 == 540


def test_get_mod2pi_angles_in_range():
    """Get the list of equivalent within a given angular range."""
    # In a 360 deg range, have only 1 angle
    a = get_mod2pi_angles_in_range((0, 360), 180, deg=True)
    assert len(a) == 1 and a[0] == 180

    # If no angles within range have no angles
    a = get_mod2pi_angles_in_range((181, 360), 180, deg=True)
    assert len(a) == 0

    # With 720 deg range, have 2 angles
    a = get_mod2pi_angles_in_range((0, 720), 180, deg=True)
    assert len(a) == 2 and a[0] == 180 and a[1] == 540

    # With 1080 deg range, have 3 angles
    a = get_mod2pi_angles_in_range((-360, 720), 180, deg=True)
    assert len(a) == 3 and a[0] == -180 and a[1] == 180 and a[2] == 540
