from __future__ import annotations

import dxtbx
from dxtbx.model.experiment_list import ExperimentList


def test_jf16M_matches_expected_hierarchy(dials_data):
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

    def recursive_test(pg1, pg2):
        assert pg1.is_similar_to(pg2)
        for c1, c2 in zip(pg1, pg2):
            recursive_test(c1, c2)

    recursive_test(d1.hierarchy(), d2.hierarchy())
