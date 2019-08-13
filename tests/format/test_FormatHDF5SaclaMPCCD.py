from __future__ import absolute_import, division, print_function

import pytest

from dxtbx.format.FormatHDF5SaclaMPCCD import FormatHDF5SaclaMPCCD
from dxtbx.model.experiment_list import ExperimentListFactory

pytest.importorskip("h5py")


@pytest.mark.xfail(
    reason="static mask isn't set correctly when reading experiment list from dictionary"
)
# https://github.com/cctbx/dxtbx/issues/70#issuecomment-520060797
def test_static_mask(dials_data):
    master_h5 = (
        dials_data("sacla_mpccd_phase3").join("MPCCD-Phase3-21528-5images.h5").strpath
    )
    assert FormatHDF5SaclaMPCCD.understand(master_h5)

    expts_from_filename = ExperimentListFactory.from_filenames([master_h5])
    expts_from_dict = ExperimentListFactory.from_dict(expts_from_filename.to_dict())
    for expts in (expts_from_filename, expts_from_dict):
        assert len(expts) == 5
        imageset = expts[0].imageset
        assert imageset.get_format_class() == FormatHDF5SaclaMPCCD
        mask = imageset.get_mask(0)
        assert len(mask) == 8
        assert [m.count(False) for m in mask] == [24296] * 8


def test_MPCCD_Phase3_21528(dials_data):
    master_h5 = (
        dials_data("sacla_mpccd_phase3").join("MPCCD-Phase3-21528-5images.h5").strpath
    )
    assert FormatHDF5SaclaMPCCD.understand(master_h5)
    expts = ExperimentListFactory.from_filenames([master_h5])
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatHDF5SaclaMPCCD

    detector = imageset.get_detector()
    beam = imageset.get_beam()

    assert len(detector) == 8
    panel = detector[0]
    assert panel.get_pixel_size() == pytest.approx((0.05, 0.05))
    assert panel.get_image_size() == (512, 1024)
    assert panel.get_trusted_range() == (-1.0, 65535.0)
    assert panel.get_fast_axis() == pytest.approx(
        (-0.0026929852392799875, -0.9999963739086762, 0.0)
    )
    assert panel.get_slow_axis() == pytest.approx(
        (-0.9999963739086762, 0.0026929852392799875, 0.0)
    )
    assert panel.get_thickness() == pytest.approx(0.3)
    assert panel.get_mu() == 0
    assert panel.get_material() == ""
    assert panel.get_origin() == pytest.approx(
        (-1.8203499755859376, 51.6243984375, -99.5)
    )
    assert panel.get_distance() == pytest.approx(99.5)

    assert imageset.get_goniometer() is None
    assert imageset.get_scan() is None

    assert beam.get_wavelength() == pytest.approx(1.2452843833238922)
    assert beam.get_direction() == (0, 0, 1)


def test_MPCCD_RECONST_MODE(dials_data, monkeypatch):
    for MPCCD_RECONST_MODE in (0, 1):
        if MPCCD_RECONST_MODE:
            monkeypatch.setenv("MPCCD_RECONST_MODE", str(MPCCD_RECONST_MODE))
        else:
            monkeypatch.delenv("MPCCD_RECONST_MODE", raising=False)

        master_h5 = (
            dials_data("sacla_mpccd_phase3")
            .join("MPCCD-Phase3-21528-5images.h5")
            .strpath
        )
        expts = ExperimentListFactory.from_filenames([master_h5])
        imageset = expts[0].imageset
        # Horrible hack to work around format_instance caching
        imageset.reader().nullify_format_instance()
        raw_data = imageset.get_raw_data(0)
        mmm = raw_data[0].as_double().as_1d().min_max_mean()
        mmm.show()
        if MPCCD_RECONST_MODE:
            assert mmm.min == -1
            assert mmm.max == 14102
            assert mmm.mean == pytest.approx(2.93752665030144)
        else:
            assert mmm.min == 0
            assert mmm.max == 3610.0
            assert mmm.mean == pytest.approx(3.397266387939453)
