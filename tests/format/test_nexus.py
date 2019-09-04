from __future__ import absolute_import, division, print_function

import numpy
import pytest

from scitbx.array_family import flex
from dxtbx.format import nexus
from dxtbx.model import Goniometer, MultiAxisGoniometer


def test_scan_factory(mocker):
    # Mock various hdf5 objects
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


def test_beam_factory(mocker):
    # Mock the hdf5 object
    obj = mocker.MagicMock()
    wavelength = mocker.MagicMock()
    wavelength.__getitem__ = lambda self, x: 1.0
    wavelength.attrs = {"units": "angstrom"}
    obj.handle = {"incident_wavelength": wavelength}

    beam = nexus.beam_factory(obj)
    assert beam.get_wavelength() == 1.0
    assert beam.get_s0() == (-0.0, -0.0, -1.0)


def test_goniometer_factory(mocker):
    # Mock the hdf5 object
    obj = mocker.MagicMock()

    # monkeypatch this function return values
    construct_axes = mocker.patch("dxtbx.format.nexus.construct_axes")
    # single-axis gonio
    construct_axes.return_value = (
        flex.vec3_double(((1.0, 0.0, 0.0),)),  # axes
        flex.double((15.0,)),  # angles
        flex.std_string(("omega",)),  # axis_names
        0,  # scan_axis
    )
    gonio = nexus.goniometer_factory(obj)
    assert isinstance(gonio, Goniometer)
    assert gonio.get_rotation_axis() == (1.0, 0.0, 0.0)

    # multi-axis gonio
    construct_axes.return_value = (
        flex.vec3_double(
            [(1.0, -0.0037, 0.002), (-0.0046, 0.0372, -0.9993), (1.0, 0.0, 0.0)]
        ),  # axes
        flex.double((0.0, 45, 90)),  # angles
        flex.std_string(("phi", "chi", "omega")),  # axis_names
        2,  # scan_axis
    )
    gonio = nexus.goniometer_factory(obj)
    assert isinstance(gonio, MultiAxisGoniometer)
    assert gonio.get_rotation_axis() == (1.0, 0.0, 0.0)
    assert list(gonio.get_names()) == ["phi", "chi", "omega"]
    assert gonio.get_fixed_rotation() == pytest.approx(
        (
            0.7071129787730329,
            0.7065597471858479,
            0.02765065835380525,
            -0.7066599864107413,
            0.7075120963132904,
            -0.007635258760788738,
            -0.02495794175606477,
            -0.01414062329050385,
            0.9995884872867717,
        )
    )
