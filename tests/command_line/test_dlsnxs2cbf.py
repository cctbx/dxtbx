import h5py
import numpy as np
import procrunner

from dxtbx.format.FormatCBFMiniEigerDLS16MSN160 import FormatCBFMiniEigerDLS16MSN160
from dxtbx.model.experiment_list import ExperimentListFactory


def test_dlsnxs2cbf(dials_data, tmp_path):
    screen = dials_data("thaumatin_eiger_screen")
    master = screen.join("Therm_6_1_master.h5")
    result = procrunner.run(
        ["dxtbx.dlsnxs2cbf", master, "junk_%04d.cbf"], working_directory=tmp_path
    )
    assert not result.returncode and not result.stderr

    output_files = ["junk_%04d.cbf" % j for j in (1, 2, 3)]
    for file in output_files:
        assert file.encode("latin-1") in result.stdout
        assert tmp_path.joinpath(file).is_file()

    expts = ExperimentListFactory.from_filenames(
        str(tmp_path / file) for file in output_files
    )
    assert all(
        imgset.get_format_class() == FormatCBFMiniEigerDLS16MSN160
        for imgset in expts.imagesets()
    )

    with h5py.File(master) as fh:
        for i, imgset in enumerate(expts.imagesets()):
            original = fh["/entry/data/data_000001"][i][()]
            sel = np.where(original < original.max())
            np.testing.assert_equal(
                fh["/entry/data/data_000001"][i][sel],
                imgset.get_raw_data(0)[0].as_numpy_array()[sel],
            )
