import pickle

import pytest

from scitbx.array_family import flex

from dxtbx.format.FormatHDF5SaclaMPCCD import FormatHDF5SaclaMPCCD
from dxtbx.format.image import ImageBool
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx_masking_ext import mask_untrusted_resolution_range

pytest.importorskip("h5py")


def test_static_mask(dials_data):
    master_h5 = (
        dials_data("image_examples").join("SACLA-MPCCD-Phase3-21528-5images.h5").strpath
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

    imageset.reader().nullify_format_instance()


def test_MPCCD_Phase3_21528(dials_data):
    master_h5 = (
        dials_data("image_examples").join("SACLA-MPCCD-Phase3-21528-5images.h5").strpath
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

    assert beam.get_wavelength() == pytest.approx(1.2452863763694901)
    assert beam.get_sample_to_source_direction() == (0, 0, 1)
    imageset.reader().nullify_format_instance()


def test_MPCCD_RECONST_MODE(dials_data, monkeypatch):
    for MPCCD_RECONST_MODE in (0, 1):
        if MPCCD_RECONST_MODE:
            monkeypatch.setenv("MPCCD_RECONST_MODE", str(MPCCD_RECONST_MODE))
        else:
            monkeypatch.delenv("MPCCD_RECONST_MODE", raising=False)

        master_h5 = (
            dials_data("image_examples")
            .join("SACLA-MPCCD-Phase3-21528-5images.h5")
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
    # Horrible hack to work around format_instance caching
    # This is needed to prevent an incorrect format_instance affecting
    # other tests that are run after this one.
    imageset.reader().nullify_format_instance()


def test_combine_with_user_static_mask(dials_data, tmpdir):
    master_h5 = (
        dials_data("image_examples").join("SACLA-MPCCD-Phase3-21528-5images.h5").strpath
    )
    assert FormatHDF5SaclaMPCCD.understand(master_h5)

    expts_from_filename = ExperimentListFactory.from_filenames([master_h5])
    for i, expt in enumerate(expts_from_filename):
        mask = []
        for panel in expt.detector:
            m = flex.bool(flex.grid(reversed(panel.get_image_size())), True)
            mask_untrusted_resolution_range(m, expt.beam, panel, 0, 2.23)
            mask.append(m)
        mask = tuple(mask)
        # exact number of masked pixels varies with wavelength between experiments
        assert [_.count(False) for _ in mask] == pytest.approx(
            [50643, 0, 0, 48003, 48356, 0, 0, 52191], abs=1e3
        )
        mask_file = tmpdir.join("pixel_%i.mask" % i)
        with mask_file.open("wb") as f:
            pickle.dump(mask, f)
        expt.imageset.external_lookup.mask.filename = mask_file.strpath
        expt.imageset.external_lookup.mask.data = ImageBool(mask)

    expts_from_dict = ExperimentListFactory.from_dict(expts_from_filename.to_dict())
    imageset = expts_from_dict[0].imageset
    mask = imageset.get_mask(0)
    assert len(mask) == 8
    assert [_.count(False) for _ in mask] == [
        65332,
        24296,
        24296,
        63335,
        63660,
        24296,
        24296,
        66691,
    ]
    imageset.reader().nullify_format_instance()


@pytest.mark.parametrize(("clear_cache"), [False, True])
def test_HDF5_format_caching(dials_data, clear_cache):
    """
    xfail: see https://github.com/cctbx/dxtbx/issues/245
    """
    img_file = "SACLA-MPCCD-Phase3-21528-5images.h5"
    master_h5 = dials_data("image_examples").join(img_file).strpath

    expts1 = ExperimentListFactory.from_filenames([master_h5])
    expts1[0].imageset.get_mask(0)
    if clear_cache:
        expts1[0].imageset.reader().nullify_format_instance()
    expts2 = ExperimentListFactory.from_filenames([master_h5])

    for e1, e2 in zip(expts1, expts2):
        assert (
            e1.imageset.get_beam().get_wavelength()
            == e2.imageset.get_beam().get_wavelength()
        )
