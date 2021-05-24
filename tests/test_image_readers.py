import os

import pycbf

from dxtbx.format.image import cbf_read_buffer
from dxtbx.model.detector import DetectorFactory


def test_cbf_buffer(dials_regression):
    filename = os.path.join(
        dials_regression, "image_examples", "dials-190", "whatev1_01_00001.cbf"
    )
    with open(filename, "rb") as f:
        contents = f.read()

    handle = pycbf.cbf_handle_struct()
    cbf_read_buffer(handle, contents, pycbf.MSG_DIGEST)
    det = DetectorFactory.imgCIF_H(handle, "unknown")
    assert det
