from __future__ import annotations

import os.path
import shutil
import subprocess


def test_plot_models_noninteractive_to_pdf(dials_data, tmp_path):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = str(ssx / "imported_no_ref_5.expt")
    result = subprocess.run(
        [shutil.which("dxtbx.plot_detector_models"), expts, "pdf_file=plot_test.pdf"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert os.path.exists(tmp_path / "plot_test.pdf")
