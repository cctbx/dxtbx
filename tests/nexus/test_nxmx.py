import h5py
import numpy as np
import pytest

from dxtbx.nexus import nxmx


def test_nxentry(nxmx_example):
    nxentry = nxmx.NXentry(nxmx_example["/entry"])
    assert nxentry.definition == "NXmx"

    assert len(nxentry.samples) == 1
    assert isinstance(nxentry.samples[0], nxmx.NXsample)
    assert len(nxentry.instruments) == 1
    assert isinstance(nxentry.instruments[0], nxmx.NXinstrument)
    assert isinstance(nxentry.source, nxmx.NXsource)
    assert len(nxentry.data) == 1
    assert isinstance(nxentry.data[0], nxmx.NXdata)


def test_nxmx(nxmx_example):
    nx = nxmx.NXmx(nxmx_example)
    assert len(nx) == 1
    assert nx.keys() == nxmx_example.keys()
    entries = nx.entries
    assert len(entries) == 1
    nxentry = entries[0]
    assert nxentry.definition == "NXmx"
    assert nxentry.path == "/entry"

    samples = nxentry.samples
    assert len(samples) == 1
    sample = samples[0]
    assert sample.name == "mysample"
    assert sample.depends_on.path == "/entry/sample/transformations/phi"
    assert sample.temperature is None
    assert sample.path == "/entry/sample"

    transformations = sample.transformations
    assert len(transformations) == 1
    axes = transformations[0].axes
    assert len(axes) == 3
    assert set(axes.keys()) == {"chi", "omega", "phi"}
    phi_depends_on = axes["phi"].depends_on
    assert phi_depends_on.path == "/entry/sample/transformations/chi"

    assert len(nxentry.instruments) == 1
    instrument = nxentry.instruments[0]
    assert instrument.name == "DIAMOND BEAMLINE I03"
    assert instrument.short_name == "I03"

    assert len(instrument.beams) == 1
    beam = instrument.beams[0]
    assert beam.incident_wavelength.to("angstrom").magnitude == 0.976223

    assert len(instrument.detectors) == 1
    detector = instrument.detectors[0]
    assert detector.description == "Eiger 16M"
    assert detector.sensor_material == "Silicon"
    assert detector.sensor_thickness.to("mm").magnitude == 0.45
    assert (
        detector.depends_on.path == "/entry/instrument/detector/transformations/det_z"
    )

    assert len(detector.modules) == 1
    module = detector.modules[0]
    assert np.all(module.data_origin == [0, 0])
    assert np.all(module.data_size == [4362, 4148])
    assert module.fast_pixel_direction.matrix.shape == (1, 4, 4)
    assert list(module.fast_pixel_direction.matrix.flatten()) == [
        1,
        0,
        0,
        -0.075,
        0,
        1,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        1,
    ]

    assert nxentry.source.name == "Diamond"
    assert nxentry.source.short_name == "DLS"


def test_get_rotation_axes(nxmx_example):
    sample = nxmx.NXmx(nxmx_example).entries[0].samples[0]
    dependency_chain = nxmx.get_dependency_chain(sample.depends_on)
    axes = nxmx.get_rotation_axes(dependency_chain)
    assert np.all(axes.is_scan_axis == [False, False, True])
    assert np.all(axes.names == ["phi", "chi", "omega"])
    assert np.all(axes.angles == [0.0, 0.0, 0.0])
    assert np.all(
        axes.axes == np.array([[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]])
    )


@pytest.mark.parametrize(
    "scan_data", [np.array(0), np.array([0])], ids=["scalar", "vector"]
)
def test_get_rotation_axis_scalar_or_vector(scan_data):
    """
    Test that single-valued rotation axis positions can be scalar or vector.

    A rotation axis with a single angular position may be recorded in a HDF5 NeXus
    file either as an array data set with a single entry, or as a scalar data set.
    Both are equally valid.  Check that they are handled correctly in get_rotation_axis.
    """
    # Create a basic h5py data set.  A non-empty string file name is required,
    # even though there is no corresponding file.
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        # Create a single data set representing the goniometer axis.
        scan_axis = f.create_dataset("dummy_axis", data=scan_data)
        # Add the attributes of a rotation scan axis aligned with the x axis.
        scan_axis.attrs["transformation_type"] = "rotation"
        scan_axis.attrs["vector"] = (1, 0, 0)
        scan_axis.attrs["units"] = "degrees"

        # Test that we can interpret the rotation axis datum.
        scan_axes = [nxmx.NXtransformationsAxis(scan_axis)]
        nxmx.get_rotation_axes(scan_axes)


def test_get_dependency_chain(nxmx_example):
    sample = nxmx.NXmx(nxmx_example).entries[0].samples[0]
    dependency_chain = nxmx.get_dependency_chain(sample.depends_on)
    assert [d.path for d in dependency_chain] == [
        "/entry/sample/transformations/phi",
        "/entry/sample/transformations/chi",
        "/entry/sample/transformations/omega",
    ]
    assert (
        str(dependency_chain)
        == """\
/entry/sample/transformations/phi = [0] degree
  @transformation_type = rotation
  @vector = [-1.  0.  0.]
  @offset = None
  @depends_on = /entry/sample/transformations/chi
/entry/sample/transformations/chi = 0 degree
  @transformation_type = rotation
  @vector = [0 0 1]
  @offset = None
  @depends_on = /entry/sample/transformations/omega
/entry/sample/transformations/omega = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9] degree
  @transformation_type = rotation
  @vector = [-1.  0.  0.]
  @offset = None
  @depends_on = ."""
    )


@pytest.fixture
def detector_depends_on_example():
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        module = f.create_group("/entry/instrument/detector/module")
        module.attrs["NX_class"] = "NXdetector_module"

        fast_pixel_direction = module.create_dataset(
            "fast_pixel_direction", data=7.5e-5
        )
        fast_pixel_direction.attrs["transformation_type"] = "translation"
        fast_pixel_direction.attrs[
            "depends_on"
        ] = "/entry/instrument/detector/module/module_offset"
        fast_pixel_direction.attrs["vector"] = np.array([-1.0, 0.0, 0.0])
        fast_pixel_direction.attrs["offset"] = np.array([0.0, 0.0, 0.0])
        fast_pixel_direction.attrs["offset_units"] = "m"
        fast_pixel_direction.attrs["units"] = "m"

        module_offset = module.create_dataset("module_offset", data=0)
        module_offset.attrs["transformation_type"] = "translation"
        module_offset.attrs[
            "depends_on"
        ] = "/entry/instrument/detector/transformations/det_z"
        module_offset.attrs["vector"] = np.array([1.0, 0.0, 0.0])
        module_offset.attrs["offset"] = np.array([0.155985, 0.166904, -0])
        module_offset.attrs["offset_units"] = "m"
        module_offset.attrs["units"] = "m"

        transformations = f.create_group("/entry/instrument/detector/transformations")
        det_z = transformations.create_dataset("det_z", data=np.array([289.3]))
        det_z.attrs["depends_on"] = b"."
        det_z.attrs["transformation_type"] = b"translation"
        det_z.attrs["units"] = b"mm"
        det_z.attrs["vector"] = np.array([0.0, 0.0, 1.0])

        yield f


def test_get_dependency_chain_detector(detector_depends_on_example):
    fast_pixel_direction = detector_depends_on_example[
        "/entry/instrument/detector/module/fast_pixel_direction"
    ]
    fast_axis = nxmx.NXtransformationsAxis(fast_pixel_direction)
    dependency_chain = nxmx.get_dependency_chain(fast_axis)
    assert len(dependency_chain) == 3
    assert [d.path for d in dependency_chain] == [
        "/entry/instrument/detector/module/fast_pixel_direction",
        "/entry/instrument/detector/module/module_offset",
        "/entry/instrument/detector/transformations/det_z",
    ]
    A = nxmx.get_cumulative_transformation(dependency_chain)
    assert A.shape == (1, 4, 4)
    assert np.allclose(
        A[0],
        np.array(
            [
                [1.0, 0.0, 0.0, 155.91],
                [0.0, 1.0, 0.0, 166.904],
                [0.0, 0.0, 1.0, 289.3],
                [0.0, 0.0, 0.0, 1.0],
            ]
        ),
    )


def test_get_cumulative_transformation(nxmx_example):
    sample = nxmx.NXmx(nxmx_example).entries[0].samples[0]
    dependency_chain = nxmx.get_dependency_chain(sample.depends_on)
    A = nxmx.get_cumulative_transformation(dependency_chain)
    assert A.shape == (10, 4, 4)
    assert np.all(
        A[0]
        == np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
    )


@pytest.fixture
def detector_group():
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        entry = f.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        entry["definition"] = "NXmx"
        instrument = entry.create_group("instrument")
        instrument.attrs["NX_class"] = "NXinstrument"

        group = instrument.create_group("detector_group")
        group.attrs["NX_class"] = "NXdetector_group"

        group.create_dataset(
            "group_names",
            data=[np.string_(n) for n in ("DET", "DTL", "DTR", "DLL", "DLR")],
            dtype="S12",
        )
        group.create_dataset("group_index", data=np.array([1, 2, 3, 4, 5]))
        group.create_dataset("group_parent", data=np.array([-1, 1, 1, 1]))
        yield f


def test_nxdetector_group(detector_group):
    group = nxmx.NXmx(detector_group).entries[0].instruments[0].detector_groups[0]
    assert list(group.group_names) == ["DET", "DTL", "DTR", "DLL", "DLR"]
    assert list(group.group_index) == [1, 2, 3, 4, 5]
    assert list(group.group_parent) == [-1, 1, 1, 1]
