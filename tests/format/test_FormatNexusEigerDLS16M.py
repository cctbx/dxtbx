from __future__ import absolute_import, division, print_function

import os

import pytest

from dxtbx.format.FormatNexusEigerDLS16M import FormatNexusEigerDLS16M
from dxtbx.model.experiment_list import ExperimentListFactory

dials = pytest.importorskip("dials")
pytest.importorskip("h5py")


@pytest.mark.parametrize(
    "master_h5",
    [
        "/dls/i04/data/2019/cm23004-1/20190109/Eiger/gw/Thaum/Thau_4/Thau_4_1_master.h5",
        "/dls/i04/data/2019/cm23004-1/20190109/Eiger/gw/Thaum/Thau_4/Thau_4_1.nxs",
    ],
)
@pytest.mark.skipif(
    not os.access("/dls/i04/data/2019/cm23004-1/20190109/Eiger", os.R_OK),
    reason="Test images not available",
)
def test_rotation_scan(master_h5):
    assert FormatNexusEigerDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames(
        [master_h5], format_kwargs={"dynamic_shadowing": True}
    )
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatNexusEigerDLS16M

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
        (-166.07661632390744, 172.5371934106162, -200.0)
    )
    assert panel.get_distance() == 200

    assert len(gonio.get_axes()) == 3
    expected_axes = ((1, 0, 0), (0, 0, -1), (1, 0, 0))
    for a1, a2 in zip(gonio.get_axes(), expected_axes):
        assert a1 == pytest.approx(a2, abs=5e-2)
    assert gonio.get_scan_axis() == 2

    assert scan.get_oscillation() == (0, 0.2)
    assert scan.get_image_range() == (1, 900)

    assert beam.get_wavelength() == pytest.approx(0.979499)
    assert beam.get_s0() == pytest.approx((0, 0, -1 / beam.get_wavelength()))


@pytest.mark.skipif(
    not os.access("/dls/i04/data/2019/cm23004-1/20190109/Eiger", os.R_OK),
    reason="Test images not available",
)
def test_grid_scan():
    master_h5 = "/dls/i04/data/2019/cm23004-1/20190109/Eiger/grid/Thaum/Thau_5/Thau_5_1_master.h5"
    assert FormatNexusEigerDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames([master_h5])
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatNexusEigerDLS16M

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

    assert scan is None

    assert beam.get_wavelength() == pytest.approx(0.979499)
    assert beam.get_s0() == pytest.approx((0, 0, -1 / beam.get_wavelength()))


@pytest.mark.xfail(
    raises=AssertionError, reason="https://github.com/cctbx/dxtbx/issues/13"
)
def test_screening(dials_data):
    master_h5 = dials_data("thaumatin_eiger_screen").join("Therm_6_1_master.h5").strpath
    assert FormatNexusEigerDLS16M.understand(master_h5)

    expts = ExperimentListFactory.from_filenames([master_h5])
    assert len(expts) == 3
    imagesets = expts[0].imageset
    assert imagesets[0].get_format_class() == FormatNexusEigerDLS16M
