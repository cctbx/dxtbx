from __future__ import absolute_import, division, print_function

import os
import pytest
import six.moves.cPickle as pickle

from dxtbx.format.FormatPYunspecifiedStill import (
    FormatPYunspecifiedStill,
    FormatPYunspecifiedStillInMemory,
)
from dxtbx.model.experiment_list import ExperimentListFactory


@pytest.mark.xfail(
    reason="static mask isn't set correctly when reading experiment list from dictionary"
)
# https://github.com/cctbx/dxtbx/issues/70#issuecomment-520060797
def test_static_mask(dials_regression):
    filename = os.path.join(
        dials_regression,
        "image_examples/LCLS_CXI/shot-s00-2011-12-02T21_07Z29.723_00569.pickle",
    )
    assert FormatPYunspecifiedStill.understand(filename)

    expts_from_filename = ExperimentListFactory.from_filenames([filename])
    expts_from_dict = ExperimentListFactory.from_dict(expts_from_filename.to_dict())
    for expts in (expts_from_filename, expts_from_dict):
        assert len(expts) == 1
        imageset = expts[0].imageset
        assert imageset.get_format_class() == FormatPYunspecifiedStill
        mask = imageset.get_mask(0)
        assert len(mask) == 1
        assert mask[0].count(False) == 867109


@pytest.mark.xfail
def test_FormatPYunspecifiedStillInMemory(dials_regression):
    filename = os.path.join(
        dials_regression,
        "image_examples/LCLS_CXI/shot-s00-2011-12-02T21_07Z29.723_00569.pickle",
    )
    assert not FormatPYunspecifiedStillInMemory.understand(filename)
    with open(filename, "rb") as f:
        d = pickle.load(f)
    assert FormatPYunspecifiedStillInMemory.understand(d)
    mem_imageset = FormatPYunspecifiedStillInMemory.get_imageset(d)  # noqa F841
