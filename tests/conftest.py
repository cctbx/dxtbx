from __future__ import annotations

import h5py
import numpy as np
import pytest


def pytest_configure():
    """Incantations to create an in-memory file in h5py."""
    pytest.h5_in_memory = {"driver": "core", "backing_store": False}


@pytest.fixture
def nxmx_example():
    with h5py.File(" ", mode="w", **pytest.h5_in_memory) as f:
        entry = f.create_group("/entry")
        entry.attrs["NX_class"] = "NXentry"
        entry["definition"] = "NXmx"
        entry["start_time"] = "2021-09-10T06:54:37Z"
        entry["end_time"] = "2021-09-10T06:55:09Z"
        entry["end_time_estimated"] = "2021-09-10T06:55:09Z"

        source = entry.create_group("source")
        source.attrs["NX_class"] = "NXsource"
        source_name = source.create_dataset("name", data="Diamond")
        source_name.attrs["short_name"] = "DLS"

        instrument = entry.create_group("instrument")
        instrument.attrs["NX_class"] = "NXinstrument"
        name = instrument.create_dataset(
            "name", data=np.string_("DIAMOND BEAMLINE I03")
        )
        name.attrs["short_name"] = "I03"

        beam = instrument.create_group("beam")
        beam.attrs["NX_class"] = "NXbeam"
        beam.create_dataset("incident_beam_size", data=np.array([3e-5, 3e-5]))
        beam["incident_beam_size"].attrs["units"] = b"m"
        beam["incident_wavelength"] = 0.976223
        beam["incident_wavelength"].attrs["units"] = b"angstrom"
        beam["total_flux"] = 1e12
        beam["total_flux"].attrs["units"] = b"Hz"

        detector = instrument.create_group("detector")
        detector.attrs["NX_class"] = "NXdetector"
        detector["beam_center_x"] = 2079.79727597266
        detector["beam_center_y"] = 2225.38773853771
        detector["count_time"] = 0.00285260857097799
        detector["depends_on"] = "/entry/instrument/detector/transformations/det_z"
        detector["description"] = "Eiger 16M"
        detector["distance"] = 0.237015940260233
        detector.create_dataset("data", data=np.zeros((100, 100)))
        detector["sensor_material"] = "Silicon"
        detector["sensor_thickness"] = 0.00045
        detector["sensor_thickness"].attrs["units"] = b"m"
        detector["x_pixel_size"] = 7.5e-05
        detector["y_pixel_size"] = 7.5e-05
        detector["underload_value"] = -1
        detector["saturation_value"] = 9266
        detector["frame_time"] = 0.1
        detector["frame_time"].attrs["units"] = "s"
        detector["bit_depth_readout"] = np.array(32)

        detector_transformations = detector.create_group("transformations")
        detector_transformations.attrs["NX_class"] = "NXtransformations"
        det_z = detector_transformations.create_dataset("det_z", data=np.array([289.3]))
        det_z.attrs["depends_on"] = b"."
        det_z.attrs["transformation_type"] = b"translation"
        det_z.attrs["units"] = b"mm"
        det_z.attrs["vector"] = np.array([0.0, 0.0, 1.0])

        module = detector.create_group("module")
        module.attrs["NX_class"] = "NXdetector_module"
        module.create_dataset("data_origin", data=np.array([0.0, 0.0]))
        module.create_dataset("data_size", data=np.array([4362, 4148]))

        fast_pixel_direction = module.create_dataset(
            "fast_pixel_direction", data=7.5e-5
        )
        fast_pixel_direction.attrs["transformation_type"] = "translation"
        fast_pixel_direction.attrs[
            "depends_on"
        ] = "/entry/instrument/detector/module/module_offset"
        fast_pixel_direction.attrs["vector"] = np.array([-1.0, 0.0, 0.0])
        fast_pixel_direction.attrs["offset"] = np.array([0.0, 0.0, 0.0])
        fast_pixel_direction.attrs["offset_units"] = b"m"
        fast_pixel_direction.attrs["units"] = b"m"

        slow_pixel_direction = module.create_dataset(
            "slow_pixel_direction", data=7.5e-5
        )
        slow_pixel_direction.attrs["transformation_type"] = "translation"
        slow_pixel_direction.attrs[
            "depends_on"
        ] = "/entry/instrument/detector/module/module_offset"
        slow_pixel_direction.attrs["vector"] = np.array([0.0, -1.0, 0.0])
        slow_pixel_direction.attrs["offset"] = np.array([0.0, 0.0, 0.0])
        slow_pixel_direction.attrs["offset_units"] = b"m"
        slow_pixel_direction.attrs["units"] = b"m"

        module_offset = module.create_dataset("module_offset", data=0)
        module_offset.attrs["transformation_type"] = "translation"
        module_offset.attrs["depends_on"] = detector["depends_on"]
        module_offset.attrs["vector"] = np.array([1.0, 0.0, 0.0])
        module_offset.attrs["offset"] = np.array([0.155985, 0.166904, -0])
        module_offset.attrs["offset_units"] = b"m"
        module_offset.attrs["units"] = b"m"

        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        sample["name"] = "mysample"
        sample["depends_on"] = b"/entry/sample/transformations/phi"
        sample["temperature"] = 273
        sample["temperature"].attrs["units"] = b"K"

        transformations = sample.create_group("transformations")
        transformations.attrs["NX_class"] = "NXtransformations"
        omega = transformations.create_dataset("omega", data=np.arange(0, 1, 0.1))
        omega.attrs["depends_on"] = b"."
        omega.attrs["transformation_type"] = b"rotation"
        omega.attrs["units"] = b"deg"
        omega.attrs["vector"] = np.array([-1.0, 0.0, 0.0])
        omega.attrs["omega_offset"] = np.array([0.0, 0.0, 0.0])

        phi = transformations.create_dataset("phi", data=np.array([0.0]))
        phi.attrs["depends_on"] = b"/entry/sample/transformations/chi"
        phi.attrs["transformation_type"] = b"rotation"
        phi.attrs["units"] = b"deg"
        phi.attrs["vector"] = np.array([-1.0, 0, 0])

        chi = transformations.create_dataset("chi", data=0.0)  # scalar dataset
        chi.attrs["depends_on"] = b"/entry/sample/transformations/omega"
        chi.attrs["transformation_type"] = b"rotation"
        chi.attrs["units"] = b"deg"
        chi.attrs["vector"] = np.array([0, 0, 1])

        data = entry.create_group("data")
        data.attrs["NX_class"] = "NXdata"

        yield f
