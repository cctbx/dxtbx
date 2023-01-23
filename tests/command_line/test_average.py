from __future__ import annotations

import shutil
import subprocess

import pytest

import dxtbx


@pytest.mark.parametrize("use_mpi", [True, False])
def test_average(dials_data, tmp_path, use_mpi):
    # Only allow MPI tests if we've got MPI capabilities
    if use_mpi:
        pytest.importorskip("mpi4py")

    data = (
        dials_data("image_examples", pathlib=True) / "SACLA-MPCCD-run266702-0-subset.h5"
    )
    if use_mpi:
        command = "mpirun"
        mpargs = "-n 2 dxtbx.image_average --mpi=True".split()
    else:
        command = shutil.which("dxtbx.image_average")
        mpargs = "-n 2".split()
    result = subprocess.run(
        [command] + mpargs + "-v -a avg.cbf -s stddev.cbf -m max.cbf".split() + [data],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    h5 = dxtbx.load(data).get_detector()
    cbf = dxtbx.load(tmp_path / "avg.cbf").get_detector()

    assert h5.is_similar_to(cbf, ignore_trusted_range=True)
    assert h5[0].get_gain() == cbf[0].get_gain()
