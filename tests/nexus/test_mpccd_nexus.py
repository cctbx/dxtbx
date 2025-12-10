from __future__ import annotations

import dxtbx


def test_mpccd_nexus_gain(dials_data):
    """
    Tests SACLA MPCCD image from CXI.DB 221
    Includes parameter data_scale_factor, which accounts for a gain of 10
    """

    try:
        h5path = dials_data("image_examples") / "SACLA-MPCCD-run197287-0-nexus.h5"
    except Exception as e:
        print(type(e), str(e))
        raise
    img = dxtbx.load(h5path)

    d = img.get_detector()

    assert d[0].get_gain() == 10
