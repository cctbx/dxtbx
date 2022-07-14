from __future__ import annotations

import os

import h5py
import numpy as np
import pytest

from scitbx.array_family import flex

from dxtbx.format.FormatNXmxDLS import FormatNXmxDLS
from dxtbx.format.FormatNXmxDLS16M import FormatNXmxDLS16M
from dxtbx.masking import SmarGonShadowMasker
from dxtbx.model.experiment_list import ExperimentListFactory

dials = pytest.importorskip("dials")
pytest.importorskip("h5py")


@pytest.mark.parametrize(
    "master_h5",
    [
        "/dls/i04/data/2021/cm28182-4/TestProteinaseK/protK3/protK3_1_master.h5",
        "/dls/i04/data/2021/cm28182-4/TestProteinaseK/protK3/protK3_1.nxs",
    ],
)
@pytest.mark.skipif(
    not os.access("/dls/i04/data/2021/cm28182-3/TestProteinaseK/protK1", os.R_OK),
    reason="Test images not available",
)
def test_rotation_scan_i04(master_h5):
    assert FormatNXmxDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames(
        [master_h5], format_kwargs={"dynamic_shadowing": True}
    )
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatNXmxDLS16M

    detector = imageset.get_detector()
    gonio = imageset.get_goniometer()
    scan = imageset.get_scan()
    beam = imageset.get_beam()

    panel = detector[0]
    assert panel.get_pixel_size() == (0.075, 0.075)
    assert panel.get_image_size() == (4148, 4362)
    assert panel.get_trusted_range() == (-1, 22726)
    assert panel.get_fast_axis() == (1, 0, 0)
    assert panel.get_slow_axis() == (0, -1, 0)
    assert panel.get_origin() == pytest.approx(
        (-158.21930645001163, 166.30914199665585, -201.83999382350086)
    )
    assert panel.get_distance() == 201.83999382350086

    assert len(gonio.get_axes()) == 3
    expected_axes = ((1, 0, 0), (0, 0, -1), (1, 0, 0))
    for a1, a2 in zip(gonio.get_axes(), expected_axes):
        assert a1 == pytest.approx(a2, abs=5e-2)
    assert gonio.get_scan_axis() == 2

    assert pytest.approx(scan.get_oscillation()) == (0, 0.1)
    assert scan.get_image_range() == (1, 3600)

    assert beam.get_wavelength() == pytest.approx(0.979499130984)
    assert beam.get_s0() == pytest.approx((0, 0, -1 / beam.get_wavelength()))


@pytest.mark.parametrize(
    "master_h5",
    [
        "/dls/i03/data/2019/cm23003-4/20190911/SmarGon/rotation_calibration4/th1_Oo_Cc_Pp_1_master.h5",
        "/dls/i03/data/2019/cm23003-4/20190911/SmarGon/rotation_calibration4/th1_Oo_Cc_Pp_1.nxs",
    ],
)
@pytest.mark.skipif(
    not os.access("/dls/i03/data/2019/cm23003-4", os.R_OK),
    reason="Test images not available",
)
def test_rotation_scan_i03_2019_run_4(master_h5):
    assert FormatNXmxDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames(
        [master_h5], format_kwargs={"dynamic_shadowing": True}
    )
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatNXmxDLS16M

    gonio = imageset.get_goniometer()
    assert list(gonio.get_angles()) == pytest.approx([45.0, 45.0, 45.0])
    assert list(gonio.get_axes().as_double()) == pytest.approx(
        [1.0, -0.0025, 0.0056, -0.006, -0.0264, -0.9996, 1.0, 0.0, 0.0]
    )
    assert list(gonio.get_names()) == ["phi", "chi", "omega"]
    assert imageset.has_dynamic_mask()
    masker = imageset.masker()
    assert isinstance(masker, SmarGonShadowMasker)
    assert masker.get_mask(imageset.get_detector(), 0)[0].count(False) == 0
    masker.get_mask(imageset.get_detector(), 50)[0].count(False) == 486717
    assert masker.get_mask(imageset.get_detector(), 100)[0].count(False) == 1092226


@pytest.mark.parametrize(
    "master_h5",
    [
        "/dls/i04/data/2020/cm26459-3/20200617/bs/lres_1_master.h5",
        "/dls/i04/data/2020/cm26459-3/20200617/bs/lres_1.nxs",
    ],
)
@pytest.mark.skipif(
    not os.access("/dls/i04/data/2020/cm26459-3/20200617/bs/", os.R_OK),
    reason="Test images not available",
)
def test_masked_i04_32bit(master_h5):
    assert FormatNXmxDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames([master_h5])
    imageset = expts[0].imageset
    assert flex.max(imageset[0][0]) != 0x7FFFFFFF


@pytest.mark.parametrize(
    "master_h5",
    [
        "/dls/i03/data/2020/cm26458-3/20200617/test_1_master.h5",
        "/dls/i03/data/2020/cm26458-3/20200617/test_1.nxs",
    ],
)
@pytest.mark.skipif(
    not os.access("/dls/i03/data/2020/cm26458-3/20200617", os.R_OK),
    reason="Test images not available",
)
def test_masked_i03_16bit(master_h5):
    assert FormatNXmxDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames([master_h5])
    imageset = expts[0].imageset
    assert flex.min(imageset[0][0]) == -1.0
    assert flex.max(imageset[0][0]) != 0xFFFF


@pytest.mark.skipif(
    not os.access("/dls/i04/data/2019/cm23004-1/20190109/Eiger", os.R_OK),
    reason="Test images not available",
)
def test_grid_scan_i04():
    master_h5 = "/dls/i04/data/2019/cm23004-1/20190109/Eiger/grid/Thaum/Thau_5/Thau_5_1_master.h5"
    assert FormatNXmxDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames([master_h5])
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatNXmxDLS16M

    detector = imageset.get_detector()
    gonio = imageset.get_goniometer()
    scan = imageset.get_scan()
    beam = imageset.get_beam()

    panel = detector[0]
    assert panel.get_pixel_size() == (0.075, 0.075)
    assert panel.get_image_size() == (4148, 4362)
    assert panel.get_trusted_range() == (-1, 65535)
    assert panel.get_fast_axis() == (1, 0, 0)
    assert panel.get_slow_axis() == (0, -1, 0)
    assert panel.get_origin() == pytest.approx(
        (-167.44717577120824, 172.46833023184868, -350.0)
    )
    assert panel.get_distance() == 350

    assert len(gonio.get_axes()) == 3
    expected_axes = ((1, 0, 0), (0, 0, -1), (1, 0, 0))
    for a1, a2 in zip(gonio.get_axes(), expected_axes):
        assert a1 == pytest.approx(a2, abs=5e-2)
    # assert gonio.get_scan_axis() == 2

    if scan:
        osc = scan.get_oscillation()
        assert osc[0] == osc[1]

    assert beam.get_wavelength() == pytest.approx(0.979499)
    assert beam.get_s0() == pytest.approx((0, 0, -1 / beam.get_wavelength()))


@pytest.mark.xfail(
    raises=AssertionError, reason="https://github.com/cctbx/dxtbx/issues/13"
)
def test_screening(dials_data):
    master_h5 = (
        dials_data("thaumatin_eiger_screen", pathlib=True) / "Therm_6_1_master.h5"
    )
    assert FormatNXmxDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames([master_h5])
    assert len(expts) == 3
    imagesets = expts[0].imageset
    assert imagesets[0].get_format_class() == FormatNXmxDLS16M


@pytest.mark.parametrize("beamline", ["I03", "I04"])
def test_understand(beamline, tmp_path):
    # See https://jira.diamond.ac.uk/browse/MXGDA-3624
    nxs = tmp_path / "data.nxs"
    with h5py.File(nxs, mode="w") as fh:
        entry = fh.create_group("entry")
        instrument = entry.create_group("instrument")
        instrument.attrs["short_name"] = np.string_(f"DLS {beamline}")
        name = instrument.create_dataset(
            "name", data=np.string_(f"DIAMOND BEAMLINE {beamline}")
        )
        name.attrs["short_name"] = np.string_(f"DLS {beamline}")
    assert FormatNXmxDLS16M.understand(nxs)
    assert FormatNXmxDLS.understand(nxs)


@pytest.mark.parametrize("beamline", ["I03", "I04"])
def test_understand_legacy(beamline, tmp_path):
    # See https://jira.diamond.ac.uk/browse/MXGDA-3624
    nxs = tmp_path / "data.nxs"
    with h5py.File(nxs, mode="w") as fh:
        entry = fh.create_group("entry")
        instrument = entry.create_group("instrument")
        instrument.attrs["short_name"] = np.string_(f"{beamline}")
        name = instrument.create_dataset("name", data=np.string_(f"{beamline}"))
        name.attrs["short_name"] = np.string_(f"{beamline}")
    assert FormatNXmxDLS16M.understand(nxs)
    assert FormatNXmxDLS.understand(nxs)


def test_do_not_understand_name_none(tmp_path):
    nxs = tmp_path / "data.nxs"
    with h5py.File(nxs, mode="w") as fh:
        entry = fh.create_group("entry")
        entry.create_group("instrument")
    assert not FormatNXmxDLS16M.understand(nxs)


def test_do_not_understand_i24(tmp_path):
    nxs = tmp_path / "data.nxs"
    with h5py.File(nxs, mode="w") as fh:
        entry = fh.create_group("entry")
        instrument = entry.create_group("instrument")
        instrument.attrs["short_name"] = np.string_("DLS I24")
        name = instrument.create_dataset(
            "name", data=np.string_("DIAMOND BEAMLINE I24")
        )
        name.attrs["short_name"] = np.string_("DLS I24")
    assert not FormatNXmxDLS16M.understand(nxs)
