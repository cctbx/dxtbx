from __future__ import absolute_import, division, print_function

import h5py
import mock
import numpy
import pytest

from scitbx.array_family import flex
from scitbx import matrix
from dxtbx.format import nexus
from dxtbx.model import BeamFactory, Goniometer, MultiAxisGoniometer


def mock_hdf5_dataset(value, **kwargs):
    dataset = mock.MagicMock(spec=h5py.Dataset, **kwargs)
    dataset.__getitem__.return_value = value
    return dataset


@pytest.fixture
def mock_find_goniometer_rotation(mocker):
    def _find_goniometer_rotation(value):
        rotation = mocker.MagicMock()
        rotation.name = "/entry/sample/transformations/omega"
        rotation.value = numpy.array(value)
        rotation.__len__.side_effect = rotation.value.__len__
        rotation.__getitem__.side_effect = rotation.value.__getitem__
        rotation.__iter__.side_effect = rotation.value.__iter__

        # monkeypatch these function return values
        find_goniometer_rotation = mocker.patch(
            "dxtbx.format.nexus.find_goniometer_rotation"
        )
        find_goniometer_rotation.return_value = rotation

    return _find_goniometer_rotation


def test_scan_factory_no_frame_time(mocker, mock_find_goniometer_rotation):
    obj = mocker.Mock()
    detector_obj = mocker.Mock()
    mock_find_goniometer_rotation((0.0, 0.1, 0.2, 0.3))

    detector_obj.handle = {}
    model = nexus.scan_factory(obj, detector_obj)
    assert len(model) == 4
    assert model.get_image_range() == (1, 4)
    assert model.get_oscillation() == (0.0, 0.1)
    assert list(model.get_exposure_times()) == [0.0, 0.0, 0.0, 0.0]
    assert list(model.get_epochs()) == [0.0, 0.0, 0.0, 0.0]


def test_scan_factory_frame_time(mocker, mock_find_goniometer_rotation):
    obj = mocker.Mock()
    detector_obj = mocker.Mock()
    mock_find_goniometer_rotation((0.0, 0.1, 0.2, 0.3))

    detector_obj.handle = {"frame_time": {(): 0.1}}
    model = nexus.scan_factory(obj, detector_obj)
    assert list(model.get_exposure_times()) == [0.1, 0.1, 0.1, 0.1]
    assert list(model.get_epochs()) == pytest.approx([0.0, 0.1, 0.2, 0.3])


def test_beam_factory(mocker):
    # Mock the hdf5 object
    obj = mocker.MagicMock()
    wavelength = mocker.MagicMock()
    wavelength.__getitem__.return_value = 1.0
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


def test_detector_factory(mocker):
    obj = mocker.MagicMock()
    beam = BeamFactory.simple(0.9802735610373182)

    handle = mocker.MagicMock()

    d = {
        "type": mock_hdf5_dataset("Pixel"),
        "saturation_value": mock_hdf5_dataset(65535),
        "sensor_thickness": mock_hdf5_dataset(0.00045, attrs={"units": "m"}),
        "sensor_material": mock_hdf5_dataset(numpy.string_("Silicon")),
        "x_pixel_size": mock_hdf5_dataset(7.5e-05),
        "y_pixel_size": mock_hdf5_dataset(7.5e-05),
        "beam_center_x": mock_hdf5_dataset(2214.3285049568203),
        "beam_center_y": mock_hdf5_dataset(2300.4972375614147),
    }
    handle.__getitem__.side_effect = d.__getitem__
    handle.__contains__.side_effect = d.__contains__
    handle.name = "/entry/instrument/detector"
    obj.handle = handle

    module = mocker.MagicMock()
    module_d = {
        "fast_pixel_direction": mock_hdf5_dataset(
            7.5e-05, attrs={"units": "m", "vector": numpy.array([-1.0, 0.0, 0.0])}
        ),
        "slow_pixel_direction": mock_hdf5_dataset(
            7.5e-05, attrs={"units": "m", "vector": numpy.array([0.0, -1.0, 0.0])}
        ),
        "module_offset": mock_hdf5_dataset(
            0.0, name="/entry/instrument/detector/module/module_offset"
        ),
        "data_size": mock_hdf5_dataset(numpy.array([4148, 4362])),
    }
    module.__getitem__.side_effect = module_d.__getitem__
    module.__contains__.side_effect = module_d.__contains__
    obj.modules[0].handle = module

    construct_vector = mocker.patch("dxtbx.format.nexus.construct_vector")
    construct_vector.return_value = matrix.col(
        (166.07463787176152, 172.53729281710608, 199.78346957333733)
    )

    detector = nexus.detector_factory(obj, beam)
    assert len(detector) == 1
    panel = detector[0]
    assert panel.get_thickness() == 0.45
    assert panel.get_pixel_size() == (0.075, 0.075)
    assert panel.get_image_size() == (4148, 4362)
    assert panel.get_fast_axis() == (1, 0, 0)
    assert panel.get_slow_axis() == (0, -1, 0)
    assert panel.get_origin() == pytest.approx(
        (-166.07463787176152, 172.53729281710608, -199.78346957333733)
    )
    assert panel.get_distance() == pytest.approx(199.78346957333733)
