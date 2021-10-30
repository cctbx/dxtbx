import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

"""
Test deserializing an experiment list that has single file indices while using check_format = True or False
"""


def test_split_single_image_imageset(dials_data, tmpdir):
    pytest.importorskip("h5py")
    sacla_file = dials_data("image_examples") / "SACLA-MPCCD-run266702-0-subset.h5"
    expts = ExperimentListFactory.from_filenames([sacla_file])
    assert len(expts) == 4
    subset = expts[2:3]
    assert len(subset) == 1
    assert list(subset[0].imageset.indices()) == [2]

    dumped_filename = tmpdir / "split.expt"
    subset.as_file(dumped_filename)

    subset = ExperimentListFactory.from_json_file(dumped_filename, check_format=True)
    assert len(subset) == 1
    assert list(subset[0].imageset.indices()) == [2]

    subset = ExperimentListFactory.from_json_file(dumped_filename, check_format=False)
    assert len(subset) == 1
    assert list(subset[0].imageset.indices()) == [2]
