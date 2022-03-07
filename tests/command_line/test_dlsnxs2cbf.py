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


def rationalise_whitespace(string: str):
    """
    Reformat argparse output in a common layout.

    argparse.parser.print_help() merges consecutive spaces and changes the linebreaks in
    argparse.parser.description, allowing breaks on a hyphen.  In order to compare
    expected and actual help output, rationalise a string to a common layout by
    replacing all instances of '-\n' with '-', and then replacing multiple consecutive
    whitespace characters with a single space.
    """
    return " ".join(string.replace("-\n", "-").split())


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


@pytest.mark.parametrize("remove_axis", ["phi", "kappa"], ids=["no phi", "no kappa"])
def test_dlsnxs2cbf_deleted_axis(dials_data, tmp_path, remove_axis):
    """Check that a master file without φ or κ axes processes OK."""
    screen = dials_data("four_circle_eiger", pathlib=True) / "03_CuHF2pyz2PF6b_P_O"
    master = "CuHF2pyz2PF6b_P_O_02.nxs"
    links = (path.name for path in screen.glob("*") if path.name != master)
    try:
        for name in links:
            (tmp_path / name).symlink_to(screen / name)
    except OSError:
        warnings.warn(
            "Copying files where unable to symlink. On Windows, Administrators "
            "or users with Developer Mode can create symlinks freely."
        )
        for name in links:
            shutil.copy(os.fspath(screen / name), os.fspath(tmp_path))
    shutil.copy(screen / master, tmp_path / master)

    with h5py.File(tmp_path / master, "r+") as f:
        del f[f"entry/sample/transformations/{remove_axis}"]
        # When removing φ, ensure that the sample's depends_on data set is redirected
        # to point to κ.  When removing κ, ensure that the depends_on attribute of φ is
        # redirected to point to ω.
        if remove_axis == "phi":
            depends_on_key = "entry/sample/depends_on"
            del f[depends_on_key]
            f.create_dataset(depends_on_key, data="/entry/sample/transformations/kappa")
        elif remove_axis == "kappa":
            phi = f["entry/sample/transformations/phi"]
            phi.attrs["depends_on"] = "/entry/sample/transformations/omega"

    make_cbf(tmp_path / master, tmp_path)


def test_dlsnxs2cbf_help(capsys):
    """Test that the expected help message appears."""
    with pytest.raises(SystemExit):
        run(["-h"])

    captured = rationalise_whitespace(capsys.readouterr().out)
    assert rationalise_whitespace(parser.description) in captured
