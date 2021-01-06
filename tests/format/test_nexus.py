import types

import h5py
import numpy as np
import pytest

from dxtbx.format import nexus


def test_mask_factory_no_mask(tmpdir):
    f = h5py.File(f"{tmpdir}.h5", "w")
    group = f.create_group("entry/instrument/detector")
    objects = [types.SimpleNamespace(handle=group)]
    assert nexus.mask_factory(objects) is None


def test_mask_factory_pixel_mask_applied(tmpdir):
    f = h5py.File(f"{tmpdir}.h5", "w")
    group = f.create_group("entry/instrument/detector")
    group["pixel_mask_applied"] = True
    objects = [types.SimpleNamespace(handle=group)]
    assert nexus.mask_factory(objects) is None


@pytest.mark.parametrize(
    "group_name",
    ["entry/instrument/detector", "entry/instrument/detector/detectorSpecific"],
)
def test_mask_factory_pixel_mask(group_name, tmpdir, mocker):
    f = h5py.File(f"{tmpdir}.h5", "w")
    group = f.create_group(group_name)
    mask = np.random.choice([True, False], size=(100, 120), p=[0.01, 0.99])
    group.create_dataset("pixel_mask", data=mask)
    mocker.patch.object(
        nexus,
        "get_detector_module_slices",
        return_value=[[slice(0, 100, 1), slice(0, 120, 1)]],
    )
    objects = [types.SimpleNamespace(handle=f["entry/instrument/detector"])]
    mask2 = nexus.mask_factory(objects)
    assert mask2 is not None
    assert (mask2[0].as_numpy_array() ^ mask).all()
