import h5py
import pytest

from dxtbx.format import nexus


def test_data_factory_nxs(dials_data):
    nxs_file = dials_data("vmxi_thaumatin") / "image_15799.nxs"
    with h5py.File(nxs_file) as fh:
        data = fh["/entry/data"]
        data_factory = nexus.DataFactory(nexus.NXdata(data))
        for i in (0, data["data"].shape[0] - 1):
            assert data_factory[i].size()
        # If we try and access an image beyond the size of the dataset it should fail.
        # See https://github.com/cctbx/dxtbx/issues/285
        with pytest.raises(IndexError):
            data_factory[data["data"].shape[0]]
