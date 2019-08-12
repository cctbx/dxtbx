from __future__ import absolute_import, division, print_function

import pytest

from dxtbx.format.FormatHDF5SaclaMPCCD import FormatHDF5SaclaMPCCD
from dxtbx.model.experiment_list import ExperimentListFactory

pytest.importorskip("h5py")


@pytest.mark.xfail(
    reason="static mask isn't set correctly when reading experiment list from dictionary"
)
# https://github.com/cctbx/dxtbx/issues/70#issuecomment-520060797
def test_static_mask(dials_data, tmpdir):
    master_h5 = (
        dials_data("sacla_mpccd_phase3").join("MPCCD-Phase3-21528-5images.h5").strpath
    )
    assert FormatHDF5SaclaMPCCD.understand(master_h5)

    expts_from_filename = ExperimentListFactory.from_filenames([master_h5])
    expts_from_dict = ExperimentListFactory.from_dict(expts_from_filename.to_dict())
    for expts in (expts_from_filename, expts_from_dict):
        assert len(expts) == 5
        imageset = expts[0].imageset
        assert imageset.get_format_class() == FormatHDF5SaclaMPCCD
        mask = imageset.get_mask(0)
        assert len(mask) == 8
        assert [m.count(False) for m in mask] == [24296] * 8
