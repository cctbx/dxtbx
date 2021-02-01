from __future__ import absolute_import, division, print_function

import procrunner

from scitbx.array_family import flex

from dxtbx.format.FormatCBFMiniEigerDLS16MSN160 import FormatCBFMiniEigerDLS16MSN160
from dxtbx.model.experiment_list import ExperimentListFactory


def test_dlsnxs2cbf(dials_data, tmpdir):
    screen = dials_data("thaumatin_eiger_screen")
    master = screen.join("Therm_6_1_master.h5")
    result = procrunner.run(
        ["dxtbx.dlsnxs2cbf", master, "junk_%04d.cbf"], working_directory=tmpdir
    )
    assert not result.returncode and not result.stderr

    expected_output = "\n".join("junk_%04d.cbf" % j for j in (1, 2, 3))

    for record in expected_output.split("\n"):
        assert record.strip().encode("latin-1") in result.stdout, record

    # check files on disk

    for j in (1, 2, 3):
        assert tmpdir.join(f"junk_{j:04d}.cbf").check()

    expts = ExperimentListFactory.from_filenames(
        [str(tmpdir.join(f"junk_{j+1:04d}.cbf")) for j in range(3)]
    )
    assert all(
        [
            imgset.get_format_class() == FormatCBFMiniEigerDLS16MSN160
            for imgset in expts.imagesets()
        ]
    )
    assert [flex.sum(imgset.get_raw_data(0)[0]) for imgset in expts.imagesets()] == [
        26990833,
        27447160,
        22271214,
    ]
