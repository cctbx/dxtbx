import dxtbx
from dxtbx.model.experiment_list import ExperimentList


def test_jf16M(tmp_path, dials_data):
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

    # expt file imported using FormatNexus, but data file is imported using
    # FormatNXmx. Panel group naming is different and hierarchy has different structure
    assert not d1 == d2

    for p1, p2 in zip(d1, d2):
        assert p1.is_similar_to(p2, ignore_trusted_range=True)
