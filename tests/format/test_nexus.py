from __future__ import absolute_import, division, print_function

import numpy
import pytest

from dxtbx.format import nexus


def test_scan_factory(mocker):
    # Mock various objects
    obj = mocker.Mock()
    detector_obj = mocker.Mock()
    rotation = mocker.MagicMock()
    rotation.name = "/entry/sample/transformations/omega"
    rotation.value = numpy.array((0.0, 0.1, 0.2, 0.3))
    rotation.__len__.return_value = len(rotation.value)

    def rotation_getitem(self, idx):
        return self.value[idx]

    rotation.__getitem__ = rotation_getitem

    # monkeypatch these function return values
    find_goniometer_rotation = mocker.patch(
        "dxtbx.format.nexus.find_goniometer_rotation"
    )
    find_goniometer_rotation.return_value = rotation

    # First test with no "frame_time" set
    detector_obj.handle = {}
    model = nexus.scan_factory(obj, detector_obj)
    assert len(model) == 4
    assert model.get_image_range() == (1, 4)
    assert model.get_oscillation() == (0.0, 0.1)
    assert list(model.get_exposure_times()) == [0.0, 0.0, 0.0, 0.0]
    assert list(model.get_epochs()) == [0.0, 0.0, 0.0, 0.0]

    # Now test with "frame_time" set
    detector_obj.handle = {"frame_time": {(): 0.1}}
    model = nexus.scan_factory(obj, detector_obj)
    assert list(model.get_exposure_times()) == [0.1, 0.1, 0.1, 0.1]
    assert list(model.get_epochs()) == pytest.approx([0.0, 0.1, 0.2, 0.3])
