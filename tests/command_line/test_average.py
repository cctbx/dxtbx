from __future__ import absolute_import, division, print_function

import os

import procrunner

import dxtbx

import pytest


@pytest.mark.parametrize("use_mpi", [True, False])
def test_average(dials_regression, tmpdir, use_mpi):
    # Only allow MPI tests if we've got MPI capabilities
    if use_mpi:
        pytest.importorskip("mpi4py")

    data = os.path.join(
        dials_regression,
        "image_examples",
        "SACLA_MPCCD_Cheetah",
        "run266702-0-subset.h5",
    )
    if use_mpi:
        command = "mpirun"
        mpargs = "-n 2 dxtbx.image_average".split()
    else:
        command = "dxtbx.image_average"
        mpargs = "-n 2".split()
    result = procrunner.run(
        [command] + mpargs + "-v -a avg.cbf -s stddev.cbf -m max.cbf".split() + [data],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    h5 = dxtbx.load(data).get_detector()
    cbf = dxtbx.load(tmpdir.join("avg.cbf")).get_detector()

    assert h5.is_similar_to(cbf)
    assert h5[0].get_gain() == cbf[0].get_gain()
