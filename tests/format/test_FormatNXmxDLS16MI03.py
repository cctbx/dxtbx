from __future__ import annotations

import os

import pytest

from dxtbx.format.FormatNXmxDLS16MI03 import FormatNXmxDLS16MI03
from dxtbx.model.experiment_list import ExperimentListFactory

dials = pytest.importorskip("dials")
pytest.importorskip("h5py")


@pytest.mark.parametrize(
    "master_h5,masked_count",
    [
        (
            "/dls/i03/data/2021/cm28170-2/TestProteinaseK/protk_4/protk_4_5.nxs",
            1052672,
        ),  # full detector mode
        (
            "/dls/i03/data/2021/cm28170-2/xraycentring/TestProteinaseK/protk_5/protk_5_15.nxs",
            526336,
        ),  # ROI mode
    ],
)
def test_masked_i03(master_h5, masked_count):
    if not os.access(master_h5, os.R_OK):
        pytest.skip("Test images not available")

    assert FormatNXmxDLS16MI03.understand(master_h5)
    expts = ExperimentListFactory.from_filenames([master_h5])
    assert (
        expts[0].detector[0].get_untrusted_rectangle_mask().count(False) == masked_count
    )
