from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.model.experiment_list import ExperimentListFactory


def test_single_panel(dials_data, tmpdir):
    filename = dials_data("x4wide") / "X4_wide_M1S4_2_0001.cbf"

    assert FormatCBFMiniPilatus.understand(filename.strpath)
    expts = ExperimentListFactory.from_filenames([filename.strpath])
    assert len(expts) == 1
    assert len(expts[0].detector) == 1


def test_multi_panel(dials_data, tmpdir):
    filename = dials_data("x4wide") / "X4_wide_M1S4_2_0001.cbf"

    assert FormatCBFMiniPilatus.understand(filename.strpath)
    expts = ExperimentListFactory.from_filenames(
        [filename.strpath], format_kwargs={"multi_panel": True}
    )
    assert len(expts) == 1
    assert len(expts[0].detector) == 60
