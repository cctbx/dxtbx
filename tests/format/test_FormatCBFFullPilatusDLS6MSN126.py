from __future__ import absolute_import, division, print_function

import os
import pytest

from dxtbx.format.FormatCBFFullPilatusDLS6MSN126 import FormatCBFFullPilatusDLS6MSN126
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.masking import SmarGonShadowMasker


def test_DLS_I03_smargon(dials_data, tmpdir):
    filename_gz = (
        dials_data("image_examples").join("DLS_I03_smargon_0001.cbf.gz").strpath
    )

    # The need to manually extract the file will be removed by
    # https://github.com/cctbx/dxtbx/pull/80
    import gzip

    filename = tmpdir.join(os.path.split(filename_gz[:-3])[-1]).strpath
    with open(filename, "wb") as f:
        with gzip.open(filename_gz, "rb") as fgz:
            f.write(fgz.read())

    assert FormatCBFFullPilatusDLS6MSN126.understand(filename)
    expts = ExperimentListFactory.from_filenames(
        [filename], format_kwargs={"dynamic_shadowing": True}
    )
    assert len(expts) == 1
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatCBFFullPilatusDLS6MSN126
    gonio = imageset.get_goniometer()
    assert list(gonio.get_angles()) == pytest.approx([45.0, 45.0, 45.0])
    assert list(gonio.get_axes().as_double()) == pytest.approx(
        [1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0]
    )
    assert list(gonio.get_names()) == ["GON_PHI", "GON_CHI", "GON_OMEGA"]
    assert imageset.has_dynamic_mask()
    masker = imageset.masker()
    assert isinstance(masker, SmarGonShadowMasker)
    assert masker.get_mask(imageset.get_detector(), 0)[0].count(False) == 0
    assert masker.get_mask(imageset.get_detector(), 100)[0].count(False) == 261588
