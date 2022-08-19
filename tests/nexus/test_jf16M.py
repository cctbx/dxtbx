import dxtbx
from dxtbx.model.experiment_list import ExperimentList
import pytest


@pytest.mark.xfail(
    raises=ValueError, reason="https://github.com/cctbx/dxtbx/issues/537"
)
def test_jf16M(tmp_path, dials_data):
    if dials_data._attempt_fetch("lysozyme_JF16M_4img") is None:
        return  # data not availible yet, remove after dials/data#389 is merged

    try:
        h5path = (
            dials_data("lysozyme_JF16M_4img", pathlib=True)
            / "lyso009a_0087.JF07T32V01_master_4img.h5"
        )
    except Exception as e:
        print(type(e), str(e))
        raise
    img = dxtbx.load(h5path)

    d1 = img.get_detector()

    expts_path = (
        dials_data("lysozyme_JF16M_4img", pathlib=True)
        / "lyso009a_0087.JF07T32V01_master_4img_imported.expt"
    )
    expts = ExperimentList.from_file(expts_path, check_format=False)

    d2 = expts[0].detector

    assert d1 == d2
