from __future__ import annotations

from dxtbx.util import ersatz_uuid4
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
