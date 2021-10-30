import pytest

from dxtbx.format.FormatCBFFullPilatusDLS6MSN126 import FormatCBFFullPilatusDLS6MSN126
from dxtbx.masking import SmarGonShadowMasker
from dxtbx.model.experiment_list import ExperimentListFactory


def test_DLS_I03_smargon(dials_data):
    filename = dials_data("image_examples").join("DLS_I03_smargon_0001.cbf.gz").strpath
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
