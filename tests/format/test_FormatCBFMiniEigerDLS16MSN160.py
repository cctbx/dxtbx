import pytest

from dxtbx.format.FormatCBFMiniEigerDLS16MSN160 import FormatCBFMiniEigerDLS16MSN160
from dxtbx.masking import SmarGonShadowMasker
from dxtbx.model.experiment_list import ExperimentListFactory


def test_dlsnxs2cbf_therm(dials_data):
    filename = dials_data("image_examples").join("dlsnxs2cbf_therm_0001.cbf.gz").strpath

    assert FormatCBFMiniEigerDLS16MSN160.understand(filename)
    expts = ExperimentListFactory.from_filenames(
        [filename], format_kwargs={"dynamic_shadowing": True}
    )
    assert len(expts) == 1
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatCBFMiniEigerDLS16MSN160
    assert imageset.has_dynamic_mask()
    gonio = imageset.get_goniometer()
    assert list(gonio.get_angles()) == pytest.approx([0.0, 0.0, 174.0])
    assert list(gonio.get_axes().as_double()) == pytest.approx(
        [1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0]
    )
    assert list(gonio.get_names()) == ["GON_PHI", "GON_CHI", "GON_OMEGA"]
    masker = imageset.masker()
    assert isinstance(masker, SmarGonShadowMasker)
    assert masker.get_mask(imageset.get_detector(), 0)[0].count(False) == 0
