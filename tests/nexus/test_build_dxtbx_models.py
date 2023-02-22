from __future__ import annotations

import h5py
import numpy as np
import nxmx
import pytest
import scipy.stats as stats

from cctbx import factor_ev_angstrom
from scitbx.array_family import flex

import dxtbx.model
import dxtbx.nexus


def test_get_dxtbx_goniometer_multi_axis(nxmx_example):
    sample = nxmx.NXmx(nxmx_example).entries[0].samples[0]
    gonio = dxtbx.nexus.get_dxtbx_goniometer(sample)
    assert isinstance(gonio, dxtbx.model.MultiAxisGoniometer)
    assert gonio.get_rotation_axis() == (1.0, 0.0, 0.0)
    assert list(gonio.get_angles()) == [0.0, 0.0, 0.0]
    assert list(gonio.get_axes()) == [
        (1.0, 0.0, 0.0),
        (0.0, 0.0, -1.0),
        (1.0, 0.0, 0.0),
    ]
    assert list(gonio.get_names()) == ["phi", "chi", "omega"]
    assert gonio.get_scan_axis() == 2


@pytest.fixture
def nxsample_single_axis():
    with h5py.File(" ", mode="w", **pytest.h5_in_memory) as f:
        entry = f.create_group("/entry")
        entry.attrs["NX_class"] = "NXentry"
        entry["definition"] = "NXmx"

        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        sample["name"] = "mysample"
        sample["depends_on"] = b"/entry/sample/transformations/omega"

        transformations = sample.create_group("transformations")
        transformations.attrs["NX_class"] = "NXtransformations"
        omega = transformations.create_dataset("omega", data=np.arange(0, 60, 0.1))
        omega.attrs["depends_on"] = b"."
        omega.attrs["transformation_type"] = b"rotation"
        omega.attrs["units"] = b"deg"
        omega.attrs["vector"] = np.array([0, 1, 0])

        yield f


def test_get_dxtbx_goniometer_single_axis(nxsample_single_axis):
    sample = nxmx.NXmx(nxsample_single_axis).entries[0].samples[0]
    gonio = dxtbx.nexus.get_dxtbx_goniometer(sample)
    assert isinstance(gonio, dxtbx.model.Goniometer)
    assert gonio.get_rotation_axis() == (0.0, 1.0, 0.0)
    assert gonio.get_fixed_rotation() == (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    assert gonio.get_setting_rotation() == (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)


@pytest.fixture
def nxsample_no_depends_on():
    with h5py.File(" ", mode="w", **pytest.h5_in_memory) as f:
        sample = f.create_group("/entry/sample")
        sample.attrs["NX_class"] = "NXsample"
        sample["depends_on"] = b"."
        yield f


def test_get_dxtbx_goniometer_sample_no_depends_on_returns_none(nxsample_no_depends_on):
    nxsample = nxmx.NXsample(nxsample_no_depends_on["/entry/sample"])
    assert dxtbx.nexus.get_dxtbx_goniometer(nxsample) is None


def test_get_dxtbx_scan_sample_no_depends_on_returns_none(
    nxsample_no_depends_on, mocker
):
    nxsample = nxmx.NXsample(nxsample_no_depends_on["/entry/sample"])
    nxdetector = mocker.Mock()
    assert dxtbx.nexus.get_dxtbx_scan(nxsample, nxdetector) is None


@pytest.fixture
def nxsample_gridscan():
    with h5py.File(" ", mode="w", **pytest.h5_in_memory) as f:
        entry = f.create_group("/entry")
        entry.attrs["NX_class"] = "NXentry"
        entry["definition"] = "NXmx"

        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        sample["name"] = "mysample"
        sample["depends_on"] = b"/entry/sample/transformations/phi"

        transformations = sample.create_group("transformations")
        transformations.attrs["NX_class"] = "NXtransformations"
        omega = transformations.create_dataset("omega", data=np.full(15, 60))
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

        chi = transformations.create_dataset("chi", data=np.array([0.0]))
        chi.attrs["depends_on"] = b"/entry/sample/transformations/sam_x"
        chi.attrs["transformation_type"] = b"rotation"
        chi.attrs["units"] = b"deg"
        chi.attrs["vector"] = np.array([0, 0, 1])

        sam_x = transformations.create_dataset("sam_x", data=np.full(15, 300))
        sam_x.attrs["depends_on"] = b"/entry/sample/transformations/sam_y"
        sam_x.attrs["equipment_component"] = b"goniometer"
        sam_x.attrs["transformation_type"] = b"translation"
        sam_x.attrs["units"] = b"mm"
        sam_x.attrs["vector"] = np.array([1.0, 0.0, 0.0])

        sam_y = transformations.create_dataset("sam_y", data=np.arange(0, 150, 10))
        sam_y.attrs["depends_on"] = b"/entry/sample/transformations/sam_z"
        sam_y.attrs["equipment_component"] = b"goniometer"
        sam_y.attrs["transformation_type"] = b"translation"
        sam_y.attrs["units"] = b"mm"
        sam_y.attrs["vector"] = np.array([0.0, 1.0, 0.0])

        sam_z = transformations.create_dataset("sam_z", data=0)
        sam_z.attrs["depends_on"] = b"/entry/sample/transformations/omega"
        sam_z.attrs["equipment_component"] = b"goniometer"
        sam_z.attrs["transformation_type"] = b"translation"
        sam_z.attrs["units"] = b"mm"
        sam_z.attrs["vector"] = np.array([0.0, 0.0, 1.0])

        instrument = entry.create_group("instrument")
        instrument.attrs["NX_class"] = "NXinstrument"

        detector = instrument.create_group("detector")
        detector.attrs["NX_class"] = "NXdetector"

        yield f


def test_get_dxtbx_goniometer_grid_scan(nxsample_gridscan):
    sample = nxmx.NXmx(nxsample_gridscan).entries[0].samples[0]
    gonio = dxtbx.nexus.get_dxtbx_goniometer(sample)
    assert isinstance(gonio, dxtbx.model.MultiAxisGoniometer)
    assert gonio.get_rotation_axis() == (1.0, 0.0, 0.0)
    assert list(gonio.get_angles()) == [0.0, 0.0, 60.0]
    assert list(gonio.get_names()) == ["phi", "chi", "omega"]
    assert gonio.get_scan_axis() == 0


def test_get_dxtbx_beam(nxmx_example):
    instrument = nxmx.NXmx(nxmx_example).entries[0].instruments[0]
    beam = dxtbx.nexus.CachedWavelengthBeamFactory(instrument.beams[0]).make_beam()
    assert isinstance(beam, dxtbx.model.Beam)
    assert beam.get_wavelength() == 0.976223
    assert beam.get_sample_to_source_direction() == (0.0, 0.0, 1.0)


def test_get_dxtbx_beam_array_length_1():
    with h5py.File(" ", mode="w", **pytest.h5_in_memory) as f:
        beam = f.create_group("/entry/instrument/beam")
        beam.attrs["NX_class"] = "NXbeam"
        beam["incident_wavelength"] = np.array([0.987])
        beam["incident_wavelength"].attrs["units"] = b"angstrom"

        nxbeam = nxmx.NXbeam(f["/entry/instrument/beam"])
        beam_factory = dxtbx.nexus.CachedWavelengthBeamFactory(nxbeam)
        assert beam_factory.make_beam(index=0).get_wavelength() == 0.987


def test_get_dxtbx_beam_array():
    wavelengths = [0.9871, 0.9872, 0.9870]
    with h5py.File(" ", mode="w", **pytest.h5_in_memory) as f:
        beam = f.create_group("/entry/instrument/beam")
        beam.attrs["NX_class"] = "NXbeam"
        beam["incident_wavelength"] = np.array(wavelengths)
        beam["incident_wavelength"].attrs["units"] = b"angstrom"

        nxbeam = nxmx.NXbeam(f["/entry/instrument/beam"])
        beam_factory = dxtbx.nexus.CachedWavelengthBeamFactory(nxbeam)
        for i, w in enumerate(wavelengths):
            assert beam_factory.make_beam(index=i).get_wavelength() == w


def test_get_dxtbx_spectrum():
    # gaussian bandpass with 20eV std dev centered on 9500 eV
    energy = 9500
    channels = np.linspace(9400, 9600, 1000)
    weights = stats.norm.pdf(channels, energy, 20)

    with h5py.File(" ", mode="w", **pytest.h5_in_memory) as f:
        beam = f.create_group("/entry/instrument/beam")
        beam.attrs["NX_class"] = "NXbeam"
        beam["incident_wavelength"] = factor_ev_angstrom / channels
        beam["incident_wavelength"].attrs["units"] = b"angstrom"
        beam["incident_wavelength_weights"] = weights

        nxbeam = nxmx.NXbeam(f["/entry/instrument/beam"])
        beam_factory = dxtbx.nexus.CachedWavelengthBeamFactory(nxbeam)
        assert (
            pytest.approx(beam_factory.make_beam().get_wavelength())
            == factor_ev_angstrom / energy
        )
        assert (
            pytest.approx(beam_factory.make_spectrum().get_weighted_energy_eV())
            == energy
        )


def test_get_dxtbx_spectrum_with_variants():
    # drifting gaussian bandpass with 20eV std dev
    energies = np.array([9500, 9520, 9510])
    channels = np.linspace(9400, 9600, 1000)
    weights = np.vstack([stats.norm.pdf(channels, e, 20) for e in energies])
    # calibrated energies not necessarily matching the weighted mean of the spectra
    calibrated_energies = np.array([e + 5 for e in energies])

    with h5py.File(" ", mode="w", **pytest.h5_in_memory) as f:
        beam = f.create_group("/entry/instrument/beam")
        beam.attrs["NX_class"] = "NXbeam"
        beam["incident_wavelength"] = factor_ev_angstrom / calibrated_energies
        beam["incident_wavelength"].attrs["units"] = b"angstrom"
        beam["incident_wavelength"].attrs["variant"] = b"incident_wavelength_1Dspectrum"
        beam["incident_wavelength_1Dspectrum"] = factor_ev_angstrom / channels
        beam["incident_wavelength_1Dspectrum"].attrs["units"] = b"angstrom"
        beam["incident_wavelength_1Dspectrum_weights"] = weights

        nxbeam = nxmx.NXbeam(f["/entry/instrument/beam"])
        beam_factory = dxtbx.nexus.CachedWavelengthBeamFactory(nxbeam)
        for i, (calibrated_energy, energy) in enumerate(
            zip(calibrated_energies, energies)
        ):
            assert (
                pytest.approx(beam_factory.make_beam(i).get_wavelength())
                == factor_ev_angstrom / calibrated_energy
            )
            assert (
                pytest.approx(beam_factory.make_spectrum(i).get_weighted_energy_eV())
                == energy
            )


def test_get_dxtbx_scan(nxmx_example):
    sample = nxmx.NXmx(nxmx_example).entries[0].samples[0]
    instrument = nxmx.NXmx(nxmx_example).entries[0].instruments[0]
    scan = dxtbx.nexus.get_dxtbx_scan(sample, instrument.detectors[0])
    assert scan.get_num_images() == 10
    assert scan.get_image_range() == (1, 10)
    assert pytest.approx(scan.get_oscillation()) == (0.0, 0.1)
    assert pytest.approx(scan.get_oscillation_range()) == (0.0, 1.0)
    assert list(scan.get_exposure_times()) == [0.1] * 10


def test_get_dxtbx_scan_grid_scan(nxsample_gridscan):
    sample = nxmx.NXmx(nxsample_gridscan).entries[0].samples[0]
    instrument = nxmx.NXmx(nxsample_gridscan).entries[0].instruments[0]
    scan = dxtbx.nexus.get_dxtbx_scan(sample, instrument.detectors[0])
    assert scan.get_num_images() == 15
    assert scan.get_image_range() == (1, 15)
    assert scan.get_oscillation() == (0.0, 0.0)
    assert scan.get_oscillation_range() == (0.0, 0.0)
    assert list(scan.get_exposure_times()) == [0.0] * 15


def test_get_dxtbx_detector(nxmx_example):
    instrument = nxmx.NXmx(nxmx_example).entries[0].instruments[0]
    wavelength = instrument.beams[0].incident_wavelength.to("angstrom").magnitude
    detector = dxtbx.nexus.get_dxtbx_detector(instrument.detectors[0], wavelength)

    assert isinstance(detector, dxtbx.model.Detector)
    assert len(detector) == 1
    panel = detector[0]
    assert panel.get_distance() == 289.3
    assert panel.get_origin() == (-155.985, 166.904, -289.3)
    assert panel.get_material() == "Si"
    assert panel.get_pixel_size() == (0.075, 0.075)
    assert panel.get_slow_axis() == (0.0, -1.0, 0.0)
    assert panel.get_fast_axis() == (1.0, 0.0, 0.0)
    assert panel.get_image_size() == (4148, 4362)
    assert panel.get_image_size_mm() == pytest.approx((311.09999999999997, 327.15))
    assert panel.get_name() == "/entry/instrument/detector/module"
    assert panel.get_normal() == (0.0, 0.0, -1.0)
    assert panel.get_trusted_range() == (0, 9266)
    assert panel.get_type() == "SENSOR_PAD"
    px_mm = panel.get_px_mm_strategy()
    assert px_mm.t0() == panel.get_thickness() == 0.45
    assert px_mm.mu() == panel.get_mu() == pytest.approx(3.9217189904637366)


def test_get_dxtbx_detector_beam_center_fallback(nxmx_example):
    nxmx_example["/entry/instrument/detector/module/module_offset"].attrs[
        "offset"
    ] = np.array((0, 0, 0))
    instrument = nxmx.NXmx(nxmx_example).entries[0].instruments[0]
    wavelength = instrument.beams[0].incident_wavelength.to("angstrom").magnitude
    detector = dxtbx.nexus.get_dxtbx_detector(instrument.detectors[0], wavelength)
    assert detector[0].get_origin() == pytest.approx(
        (-155.985, 166.904, -289.3), rel=1e-5
    )


@pytest.fixture
def detector_with_two_theta():
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        beam = f.create_group("/entry/instrument/beam")
        beam.attrs["NX_class"] = "NXbeam"
        beam["incident_wavelength"] = 0.495937
        beam["incident_wavelength"].attrs["units"] = b"angstrom"

        detector = f.create_group("/entry/instrument/detector")
        detector.attrs["NX_class"] = "NXdetector"
        detector["sensor_material"] = "Silicon"
        detector["sensor_thickness"] = 0.00045
        detector["sensor_thickness"].attrs["units"] = b"m"

        module = detector.create_group("module")
        module.attrs["NX_class"] = "NXdetector_module"
        module.create_dataset("data_size", data=np.array([2162, 2068]))

        fast_pixel_direction = module.create_dataset(
            "fast_pixel_direction", data=7.5e-5
        )
        fast_pixel_direction.attrs["transformation_type"] = "translation"
        fast_pixel_direction.attrs[
            "depends_on"
        ] = "/entry/instrument/detector/module/module_offset"
        fast_pixel_direction.attrs["vector"] = np.array([-1.0, 0.0, 0.0])
        fast_pixel_direction.attrs["units"] = "m"

        slow_pixel_direction = module.create_dataset(
            "slow_pixel_direction", data=7.5e-5
        )
        slow_pixel_direction.attrs["transformation_type"] = "translation"
        slow_pixel_direction.attrs[
            "depends_on"
        ] = "/entry/instrument/detector/module/module_offset"
        slow_pixel_direction.attrs["vector"] = np.array([0, -1.0, 0.0])
        slow_pixel_direction.attrs["units"] = "m"

        module_offset = module.create_dataset("module_offset", data=112.19250476301882)
        module_offset.attrs["transformation_type"] = "translation"
        module_offset.attrs[
            "depends_on"
        ] = "/entry/instrument/detector/transformations/det_z"
        module_offset.attrs["vector"] = np.array([0.72264186, 0.69122265, 0.0])
        module_offset.attrs["units"] = b"mm"

        transformations = detector.create_group("transformations")
        det_z = transformations.create_dataset("det_z", data=120.0)
        det_z.attrs[
            "depends_on"
        ] = b"/entry/instrument/detector/transformations/two_theta"
        det_z.attrs["transformation_type"] = b"translation"
        det_z.attrs["units"] = b"mm"
        det_z.attrs["vector"] = np.array([0.0, 0.0, 1.0])

        two_theta = transformations.create_dataset("two_theta", data=45)
        two_theta.attrs["depends_on"] = b"."
        two_theta.attrs["transformation_type"] = b"rotation"
        two_theta.attrs["units"] = b"deg"
        two_theta.attrs["vector"] = np.array([-1.0, 0.0, 0.0])

        yield f


def test_get_dxtbx_detector_with_two_theta(detector_with_two_theta):
    det = nxmx.NXdetector(detector_with_two_theta["/entry/instrument/detector"])
    wavelength = detector_with_two_theta["/entry/instrument/beam/incident_wavelength"][
        ()
    ]

    detector = dxtbx.nexus.get_dxtbx_detector(det, wavelength)
    panel = detector[0]
    assert panel.get_fast_axis() == (1.0, 0.0, 0.0)
    assert panel.get_slow_axis() == (0.0, -0.7071067811865475, -0.7071067811865476)
    assert panel.get_origin() == (
        -81.07500032000678,
        139.68894494331983,
        -30.01668254145155,
    )
    assert panel.get_distance() == pytest.approx(120)


@pytest.fixture(
    params=[True, False], ids=["equipment_component", "no equipment_component"]
)
def hierarchical_detector(request):
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        beam = f.create_group("/entry/instrument/beam")
        beam.attrs["NX_class"] = "NXbeam"
        beam["incident_wavelength"] = 0.495937
        beam["incident_wavelength"].attrs["units"] = b"angstrom"

        detector = f.create_group("/entry/instrument/ELE_D0")
        detector.attrs["NX_class"] = "NXdetector"
        detector["sensor_material"] = "Si"
        detector["sensor_thickness"] = 320
        detector["sensor_thickness"].attrs["units"] = b"microns"

        def make_module(name, depends_on, data_origin, fast_direction, slow_direction):
            module = detector.create_group(name)
            module.attrs["NX_class"] = "NXdetector_module"
            module.create_dataset("data_size", data=np.array([254, 254]))
            module.create_dataset("data_origin", data=np.array(data_origin))
            fast = module.create_dataset("fast_pixel_direction", data=0.075)
            fast.attrs["transformation_type"] = "translation"
            fast.attrs["depends_on"] = depends_on
            fast.attrs["vector"] = np.array(fast_direction)
            fast.attrs["units"] = "mm"
            slow = module.create_dataset("slow_pixel_direction", data=0.075)
            slow.attrs["transformation_type"] = "translation"
            slow.attrs["depends_on"] = depends_on
            slow.attrs["vector"] = np.array(slow_direction)
            slow.attrs["units"] = "mm"

        make_module(
            name="ARRAY_D0Q0M0A0",
            depends_on="/entry/instrument/ELE_D0/transformations/AXIS_D0Q0M0A0",
            data_origin=[1, 1],
            fast_direction=[-0.999998, -0.001781, 0],
            slow_direction=[-0.001781, 0.999998, 0],
        )
        make_module(
            name="ARRAY_D0Q0M0A1",
            depends_on="/entry/instrument/ELE_D0/transformations/AXIS_D0Q0M0A1",
            data_origin=[1, 259],
            fast_direction=[-0.999998, -0.001781, 0],
            slow_direction=[-0.001781, 0.999998, 0],
        )
        make_module(
            name="ARRAY_D0Q0M12A0",
            depends_on="/entry/instrument/ELE_D0/transformations/AXIS_D0Q0M12A0",
            data_origin=[6169, 1],
            fast_direction=[-0.999996, 0.00279501, 0],
            slow_direction=[0.00279501, 0.999996, 0],
        )
        make_module(
            name="ARRAY_D0Q0M12A1",
            depends_on="/entry/instrument/ELE_D0/transformations/AXIS_D0Q0M12A1",
            data_origin=[6169, 259],
            fast_direction=[-0.999996, 0.00279501, 0],
            slow_direction=[0.00279501, 0.999996, 0],
        )

        transformations = detector.create_group("transformations")

        t = transformations.create_dataset("AXIS_D0Q0M0A0", data=0)
        t.attrs["depends_on"] = b"AXIS_D0Q0M0"
        t.attrs["transformation_type"] = b"rotation"
        t.attrs["offset"] = np.array([38.5842, -19.1314, 0])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, -1])
        t.attrs["units"] = b"degrees"

        t = transformations.create_dataset("AXIS_D0Q0M0A1", data=0)
        t.attrs["depends_on"] = b"AXIS_D0Q0M0"
        t.attrs["transformation_type"] = b"rotation"
        t.attrs["offset"] = np.array([19.2342, -19.1658, 0])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, -1])
        t.attrs["units"] = b"degrees"

        t = transformations.create_dataset("AXIS_D0Q0M12A0", data=0)
        t.attrs["depends_on"] = b"AXIS_D0Q0M12"
        t.attrs["transformation_type"] = b"rotation"
        t.attrs["offset"] = np.array([38.4963, -19.3077, 0])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, -1])
        t.attrs["units"] = b"degrees"

        t = transformations.create_dataset("AXIS_D0Q0M12A1", data=0)
        t.attrs["depends_on"] = b"AXIS_D0Q0M12"
        t.attrs["transformation_type"] = b"rotation"
        t.attrs["offset"] = np.array([19.1463, -19.2536, 0])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, -1])
        t.attrs["units"] = b"degrees"

        t = transformations.create_dataset("AXIS_D0Q0M0", data=0)
        t.attrs["depends_on"] = b"AXIS_D0Q0"
        t.attrs["transformation_type"] = b"rotation"
        t.attrs["offset"] = np.array([38.3107, -61.1138, 0])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, -1])
        t.attrs["units"] = b"degrees"

        t = transformations.create_dataset("AXIS_D0Q0M12", data=0)
        t.attrs["depends_on"] = b"AXIS_D0Q0"
        t.attrs["transformation_type"] = b"rotation"
        t.attrs["offset"] = np.array([39.1843, 61.4038, 0])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, -1])
        t.attrs["units"] = b"degrees"

        t = transformations.create_dataset("AXIS_D0Q0", data=0)
        t.attrs["depends_on"] = b"AXIS_D0"
        t.attrs["transformation_type"] = b"rotation"
        t.attrs["offset"] = np.array([75.1333, -84.6597, 0])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, -1])
        t.attrs["units"] = b"degrees"

        t = transformations.create_dataset("AXIS_D0", data=0)
        t.attrs["depends_on"] = b"AXIS_RAIL"
        t.attrs["transformation_type"] = b"rotation"
        t.attrs["offset"] = np.array([-0.429033, 0.121021, -0.294])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, -1])
        t.attrs["units"] = b"degrees"
        if request.param:
            t.attrs["equipment_component"] = "detector_arm"

        t = transformations.create_dataset("AXIS_RAIL", data=97.83)
        t.attrs["depends_on"] = b"."
        t.attrs["transformation_type"] = b"translation"
        t.attrs["offset"] = np.array([0, 0, 0])
        t.attrs["offset_units"] = b"mm"
        t.attrs["vector"] = np.array([0, 0, 1])
        t.attrs["units"] = b"mm"
        if request.param:
            t.attrs["equipment_component"] = "detector_arm"

        yield f


def test_get_dxtbx_detector_hierarchical(hierarchical_detector):
    det = nxmx.NXdetector(hierarchical_detector["/entry/instrument/ELE_D0"])
    wavelength = hierarchical_detector["/entry/instrument/beam/incident_wavelength"][()]

    detector = dxtbx.nexus.get_dxtbx_detector(det, wavelength)
    assert len(detector) == 4

    assert detector[0].get_name() == "/entry/instrument/ELE_D0/ARRAY_D0Q0M0A0"
    assert detector[1].get_name() == "/entry/instrument/ELE_D0/ARRAY_D0Q0M0A1"
    assert detector[2].get_name() == "/entry/instrument/ELE_D0/ARRAY_D0Q0M12A0"
    assert detector[3].get_name() == "/entry/instrument/ELE_D0/ARRAY_D0Q0M12A1"

    assert (
        detector[0].parent().parent().get_name()
        == "/entry/instrument/ELE_D0/transformations/AXIS_D0Q0M0"
    )
    assert (
        detector[1].parent().parent().get_name()
        == "/entry/instrument/ELE_D0/transformations/AXIS_D0Q0M0"
    )
    assert (
        detector[2].parent().parent().get_name()
        == "/entry/instrument/ELE_D0/transformations/AXIS_D0Q0M12"
    )
    assert (
        detector[3].parent().parent().get_name()
        == "/entry/instrument/ELE_D0/transformations/AXIS_D0Q0M12"
    )
    for i in range(4):
        assert (
            detector[i].parent().parent().parent().get_name()
            == "/entry/instrument/ELE_D0/transformations/AXIS_D0Q0"
        )
        assert (
            detector[i].parent().parent().parent().parent().get_name()
            == "/entry/instrument/ELE_D0/transformations/AXIS_D0"
        )

    assert detector[0].get_origin() == (-151.599167, -164.783879, -97.536)
    assert detector[0].get_local_origin() == (0.0, 0.0, 0.0)
    fast_axis = pytest.approx((0.9999984140169291, -0.0017810007373656254, 0.0))
    assert detector[0].get_fast_axis() == fast_axis
    assert detector[0].get_local_fast_axis() == fast_axis
    slow_axis = pytest.approx((0.0017810007373656254, 0.9999984140169291, 0.0))
    assert detector[0].get_slow_axis() == slow_axis
    assert detector[0].get_local_slow_axis() == slow_axis
    assert detector[0].get_distance() == pytest.approx(97.536)

    depth = 0
    pg = detector.hierarchy()
    while True:
        depth += 1
        if pg.is_panel():
            break
        else:
            pg = pg[0]
    if (
        "equipment_component"
        in hierarchical_detector[
            "/entry/instrument/ELE_D0/transformations/AXIS_RAIL"
        ].attrs
    ):
        assert depth == 5
    else:
        assert depth == 6


@pytest.fixture
def pixel_mask_example():
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        detector = f.create_group("/entry/instrument/detector")
        detector.attrs["NX_class"] = "NXdetector"

        module = detector.create_group("module")
        module.attrs["NX_class"] = "NXdetector_module"
        module.create_dataset("data_origin", data=np.array([0.0, 0.0]))
        module.create_dataset("data_size", data=np.array([4362, 4148]))

        nf = 1028  # module pixels fast
        ns = 512  # module pixels slow
        df = 12  # module gap fast
        ds = 38  # module gap slow

        mask = np.zeros((4362, 4148), dtype="i8")
        for j in range(mask.shape[1] // (nf + df)):
            mask[:, (j + 1) * nf + j * df : (j + 1) * (nf + df)] = 1
        for i in range(mask.shape[0] // (ns + ds)):
            mask[(i + 1) * ns + i * ds : (i + 1) * (ns + ds), :] = 1

        detector.create_dataset("pixel_mask", data=mask)
        yield f


def test_get_static_mask(pixel_mask_example):
    det = nxmx.NXdetector(pixel_mask_example["/entry/instrument/detector"])
    mask = dxtbx.nexus.get_static_mask(det)
    assert len(mask) == 1
    assert isinstance(mask[0], flex.bool)
    assert mask[0].all() == (4362, 4148)
    assert np.all(
        mask[0].as_numpy_array()
        == (pixel_mask_example["/entry/instrument/detector/pixel_mask"][()] == 0)
    )


@pytest.fixture
def nxdata_example():
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        detector = f.create_group("/entry/instrument/detector")
        detector.attrs["NX_class"] = "NXdetector"

        module = detector.create_group("module")
        module.attrs["NX_class"] = "NXdetector_module"
        module.create_dataset("data_origin", data=np.array([0.0, 0.0]))
        module.create_dataset("data_size", data=np.array([4362, 4148]))

        nxdata = f.create_group("/entry/data")
        nxdata.attrs["NX_class"] = "NXdata"
        nxdata.create_dataset(
            "data",
            data=np.array([np.full((4362, 4148), i, dtype=np.int32) for i in range(3)]),
        )
        nxdata.attrs["signal"] = "/entry/data/data"

        yield f


def test_get_raw_data_single_panel(nxdata_example):
    det = nxmx.NXdetector(nxdata_example["/entry/instrument/detector"])
    nxdata = nxmx.NXdata(nxdata_example["/entry/data"])
    for i in range(3):
        raw_data = dxtbx.nexus.get_raw_data(nxdata, det, i)
        assert len(raw_data) == 1
        assert isinstance(raw_data[0], flex.int)
        assert raw_data[0].all() == (4362, 4148)
        assert raw_data[0].all_eq(i)


def test_dataset_as_flex_int():
    slices = ()
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        g = f.create_group("/foo")
        g.create_dataset("int32", data=np.array([0, 1], dtype=np.int32))
        g.create_dataset("intc", data=np.array([1, 2], dtype=np.intc))
        g.create_dataset("int64", data=np.array([2, 3], dtype=np.int64))
        for k in ("int32", "intc"):
            flex_a = dxtbx.nexus._dataset_as_flex(g[k], slices)
            assert isinstance(flex_a, flex.int)
            assert flex_a.all() == g[k].shape
            assert list(flex_a) == list(g[k])
        with pytest.raises(TypeError, match="Unsupported integer dtype .*"):
            dxtbx.nexus._dataset_as_flex(g["int64"], slices)
        flex_a = dxtbx.nexus._dataset_as_flex(g["int64"], slices, bit_depth=32)
        assert isinstance(flex_a, flex.int)
        assert flex_a.all() == g["int64"].shape
        assert list(flex_a) == list(g["int64"])


def test_dataset_as_flex_float():
    slices = ()
    np_float_types = (np.half, np.single, np.float16, np.float32)
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        g = f.create_group("/foo")
        for i, dtype in enumerate(np_float_types):
            d = g.create_dataset(f"float-{i}", data=np.array([0, 1], dtype=dtype))
            flex_a = dxtbx.nexus._dataset_as_flex(d, slices)
            assert isinstance(flex_a, flex.float)
            assert flex_a.all() == d.shape
            assert list(flex_a) == list(d)


def test_dataset_as_flex_double():
    slices = ()
    np_double_types = (
        np.float_,
        np.double,
        np.float64,
    )
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        g = f.create_group("/foo")
        for i, dtype in enumerate(np_double_types):
            d = g.create_dataset(f"double-{i}", data=np.array([0, 1], dtype=dtype))
            flex_a = dxtbx.nexus._dataset_as_flex(d, slices)
            assert isinstance(flex_a, flex.double)
            assert flex_a.all() == d.shape
            assert list(flex_a) == list(d)


def test_dataset_as_flex_unsupported():
    slices = ()
    with h5py.File(" ", "w", **pytest.h5_in_memory) as f:
        g = f.create_group("/foo")
        d = g.create_dataset("bool", data=np.array([0, 1], dtype=bool))
        with pytest.raises(TypeError, match="Unsupported dtype .*"):
            dxtbx.nexus._dataset_as_flex(d, slices)
