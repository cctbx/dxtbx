from __future__ import annotations

import os
import shutil
import warnings

import h5py
import numpy as np
import pytest

from dxtbx.command_line.dlsnxs2cbf import parser, run
from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.util.dlsnxs2cbf import make_cbf


def test_dlsnxs2cbf(dials_data, tmp_path, capsys):
    """Test basic behaviour of dxtbx.dlsnxs2cbf."""
    screen = dials_data("four_circle_eiger", pathlib=True)
    master = screen / "03_CuHF2pyz2PF6b_P_O" / "CuHF2pyz2PF6b_P_O_02.nxs"
    run([str(master), "-o", str(tmp_path)])

    output_files = [tmp_path / f"{master.stem}_{j:03d}.cbf" for j in (1, 2, 3)]
    captured = capsys.readouterr()
    assert f"{master.stem}_###.cbf" in captured.out
    assert all(path.is_file() for path in output_files)

    expts = ExperimentListFactory.from_filenames(output_files)
    assert all(
        imgset.get_format_class() == FormatCBFMiniEiger for imgset in expts.imagesets()
    )

    with h5py.File(master) as fh:
        for i, imgset in enumerate(expts.imagesets()):
            original = fh["/entry/data/data_000001"][i][()]
            sel = np.where(original < original.max())
            np.testing.assert_equal(
                original[sel], imgset.get_raw_data(0)[0].as_numpy_array()[sel]
            )


@pytest.mark.xfail(reason="The dependency chain is broken")
@pytest.mark.parametrize("remove_axis", ["phi", "chi"], ids=["no phi", "no chi"])
def test_dlsnxs2cbf_deleted_axis(dials_data, tmp_path, remove_axis):
    """Check that a master file without φ or χ axes processes OK."""
    screen = dials_data("thaumatin_eiger_screen", pathlib=True)
    master = "Therm_6_1_master.h5"
    links = (path.name for path in screen.glob("*") if path.name != master)
    try:
        for name in links:
            (tmp_path / name).symlink_to(screen / name)
    except OSError:
        warnings.warn(
            "Copying files where unable to symlink. On Windows, Administrators"
            " or users with Developer Mode can create symlinks freely."
        )
        for name in links:
            shutil.copy(os.fspath(screen / name), os.fspath(tmp_path))
    shutil.copy(screen / master, tmp_path / master)

    with h5py.File(tmp_path / master, "r+") as f:
        del f[f"entry/sample/transformations/{remove_axis}"]

    make_cbf(tmp_path / master, template=str(tmp_path / "image_%04d.cbf"))


@pytest.mark.xfail(reason="Broken for old data while collecting new data")
def test_dlsnxs2cbf_help(capsys):
    with pytest.raises(SystemExit):
        run(["-h"])
    captured = capsys.readouterr()
    assert parser.description in captured.out
    assert "Template cbf output name e.g. 'image_%04d.cbf'" in captured.out
