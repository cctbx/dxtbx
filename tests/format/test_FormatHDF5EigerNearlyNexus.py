from __future__ import absolute_import, division, print_function

import os

import pytest

from dxtbx.format.FormatHDF5EigerNearlyNexus import FormatHDF5EigerNearlyNexus
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model.goniometer import Goniometer

pytest.importorskip("h5py")


def test_semi_synthetic_dectris_eiger_nearly_nexus(dials_data, tmpdir):
    master_h5 = dials_data("image_examples").join("dectris_eiger_master.h5").strpath

    if not os.access(master_h5, os.R_OK):
        pytest.skip("Test images not available")

    assert FormatHDF5EigerNearlyNexus.understand(master_h5)

    expts = ExperimentListFactory.from_filenames([master_h5])
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatHDF5EigerNearlyNexus

    detector = imageset.get_detector()
    gonio = imageset.get_goniometer()
    scan = imageset.get_scan()
    beam = imageset.get_beam()

    panel = detector[0]
    assert panel.get_pixel_size() == pytest.approx((0.075, 0.075))
    assert panel.get_image_size() == (3110, 3269)
    assert panel.get_trusted_range() == (-1, 12440)
    assert panel.get_fast_axis() == (1, 0, 0)
    assert panel.get_slow_axis() == (0, -1, 0)
    assert panel.get_thickness() == pytest.approx(0.45)
    assert panel.get_mu() == pytest.approx(3.96763)
    assert panel.get_material() == "Si"
    assert panel.get_origin() == pytest.approx((-120.556, 118.982, -134.255), abs=1e-3)
    assert panel.get_distance() == pytest.approx(134.255)

    assert isinstance(gonio, Goniometer)
    assert gonio.get_rotation_axis() == (1, 0, 0)
    assert gonio.get_fixed_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)
    assert gonio.get_setting_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)

    assert scan.get_oscillation() == pytest.approx((0, 0))
    assert scan.get_image_range() == (1, 1)

    assert beam.get_wavelength() == pytest.approx(0.980112, abs=1e-5)
    assert beam.get_s0() == pytest.approx((0, 0, -1 / beam.get_wavelength()))


@pytest.mark.skip(reason="test needs a smaller dataset in dials-data")
def test_soleil_Proxima2A_zenodo_1443110_data03():
    # https://zenodo.org/record/1221344#.XEHr_5ynx2Q
    master_h5 = "/dls/mx-scratch/rjgildea/zenodo/1221344/200Hz/3_5_200Hz_1_master.h5"

    if not os.access(master_h5, os.R_OK):
        pytest.skip("Test images not available")

    assert FormatHDF5EigerNearlyNexus.understand(master_h5)

    expts = ExperimentListFactory.from_filenames([master_h5])
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatHDF5EigerNearlyNexus

    detector = imageset.get_detector()
    gonio = imageset.get_goniometer()
    scan = imageset.get_scan()
    beam = imageset.get_beam()

    panel = detector[0]
    assert panel.get_pixel_size() == pytest.approx((0.075, 0.075))
    assert panel.get_image_size() == (3110, 3269)
    assert panel.get_trusted_range() == (-1, 12440)
    assert panel.get_fast_axis() == (1, 0, 0)
    assert panel.get_slow_axis() == (0, -1, 0)
    assert panel.get_thickness() == pytest.approx(0.45)
    assert panel.get_mu() == pytest.approx(3.96763)
    assert panel.get_material() == "Si"
    assert panel.get_origin() == pytest.approx((-120.556, 118.982, -134.255), abs=1e-3)
    assert panel.get_distance() == pytest.approx(134.255)

    assert isinstance(gonio, Goniometer)
    assert gonio.get_rotation_axis() == (1, 0, 0)
    assert gonio.get_fixed_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)
    assert gonio.get_setting_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)

    assert scan.get_oscillation() == pytest.approx((0, 0.5))
    assert scan.get_image_range() == (1, 800)

    assert beam.get_wavelength() == pytest.approx(0.980112, abs=1e-5)
    assert beam.get_s0() == pytest.approx((0, 0, -1 / beam.get_wavelength()))
