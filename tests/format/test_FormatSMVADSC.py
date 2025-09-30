from __future__ import annotations

import bz2

import dxtbx


def test_noninteger_pedestal(dials_data, tmpdir):
    data_dir = dials_data("image_examples", pathlib=True)
    filename = str(data_dir / "APS_14BMC-q315_1_001.img.bz2")
    # Read this file in as data
    with bz2.open(filename, "rb") as f:
        data = f.read()

    # Write out with an inserted header item of a noninteger pedestal
    test_file = tmpdir / "test_pedestal_001.img"
    with test_file.open("wb") as f:
        f.write(data.replace(b"DIM=2;\n", b"DIM=2;\nIMAGE_PEDESTAL=42.0;\n"))

    # Make sure this loads
    dxtbx.load(test_file)
