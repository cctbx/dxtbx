from dxtbx.util import format_float_with_standard_uncertainty as ffwsu


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
