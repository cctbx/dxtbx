from __future__ import absolute_import, division, print_function

import os

import pycbf

from dxtbx.model.detector import DetectorFactory

try:
    from dxtbx.format.image import cbf_read_buffer
except ImportError:
    # If we didn't compile cbf_read_buffer, use the pycbf version
    def cbf_read_buffer(handle, contents, flags):
        """Use cbf_read_buffer via pycbf"""
        assert hasattr(pycbf.cbf_handle_struct, "read_buffer")
        return handle.read_buffer(contents, flags)


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
