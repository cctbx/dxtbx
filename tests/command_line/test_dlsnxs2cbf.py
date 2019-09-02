from __future__ import absolute_import, division, print_function

import os
import procrunner
import pytest  # noqa


def test_dlsnxs2cbf(dials_data, tmpdir):
    screen = dials_data("i04_eiger_screen")
    master = os.path.join(screen, "Therm_6_1_master.h5")
    result = procrunner.run(
        ["dxtbx.dlsnxs2cbf", master, "junk_%04d.cbf"], working_directory=tmpdir
    )
    assert not result.returncode and not result.stderr

    # allow extra lines to have been added (these may be comments)

    expected_output = "\n".join(["junk_%04d.cbf" % j for j in (1, 2, 3)])

    for record in expected_output.split("\n"):
        assert record.strip().encode("latin-1") in result.stdout, record
