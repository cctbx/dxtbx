from __future__ import annotations

from pathlib import Path

from dxtbx.model.beam import BeamFactory


def test_beam():
    image = str(Path(__file__).parent / "phi_scan_001.cbf")
    assert BeamFactory.imgCIF(image)
