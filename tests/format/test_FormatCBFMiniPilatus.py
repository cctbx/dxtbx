import pytest

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.model.experiment_list import ExperimentListFactory


def test_single_panel(dials_data, tmpdir):
    filename = dials_data("x4wide") / "X4_wide_M1S4_2_0001.cbf"

    assert FormatCBFMiniPilatus.understand(filename.strpath)
    expts = ExperimentListFactory.from_filenames([filename.strpath])
    assert len(expts) == 1
    assert len(expts[0].detector) == 1
    assert len(expts[0].imageset.get_raw_data(0)) == 1


def test_multi_panel(dials_data, tmpdir):
    filename = dials_data("x4wide") / "X4_wide_M1S4_2_0001.cbf"

    assert FormatCBFMiniPilatus.understand(filename.strpath)
    expts = ExperimentListFactory.from_filenames(
        [filename.strpath], format_kwargs={"multi_panel": True}
    )
    assert len(expts) == 1
    assert len(expts[0].detector) == 60
    assert len(expts[0].imageset.get_raw_data(0)) == 60
    expected_xoffset = [-212.47848, -127.51048, -42.54248, 42.42552, 127.39352]
    expected_yoffset = [
        220.00176,
        183.53776,
        147.07376,
        110.60976,
        74.14576,
        37.68176,
        1.21776,
        -35.24624,
        -71.71024,
        -108.17424,
        -144.63824,
        -181.10224,
    ]
    for i, p in enumerate(expts[0].detector):
        assert p.get_image_size() == (487, 195)
        origin = p.get_origin()
        assert origin[0] == pytest.approx(expected_xoffset[i % 5])
        assert origin[1] == pytest.approx(expected_yoffset[i // 5])
        assert origin[2] == pytest.approx(-190.18)
