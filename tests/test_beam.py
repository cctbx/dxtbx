import os

import dxtbx
from dxtbx.model.beam import BeamFactory


def test_beam():
    dxtbx_dir = dxtbx.__path__[0]

    image = os.path.join(dxtbx_dir, "tests", "phi_scan_001.cbf")
    assert BeamFactory.imgCIF(image)
