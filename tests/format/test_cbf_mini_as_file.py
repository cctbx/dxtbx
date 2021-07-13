import os

import pytest

from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.model.compare import sequence_diff
from dxtbx.model.experiment_list import (
    BeamComparison,
    DetectorComparison,
    ExperimentListFactory,
    GoniometerComparison,
)


@pytest.mark.parametrize(
    "image_file",
    [
        "image_examples/ALS_831/q315r_lyso_001.img",
        "image_examples/DLS_I02/X4_wide_M1S4_1_0001.cbf",
    ],
)
def test_cbf_writer(image_file, dials_regression, run_in_tmpdir):
    filename = os.path.join(dials_regression, image_file)
    imageset = ExperimentListFactory.from_filenames([filename])[0].imageset

    FormatCBFMini.as_file(
        imageset.get_detector(),
        imageset.get_beam(),
        imageset.get_goniometer(),
        imageset.get_sequence(),
        imageset.get_raw_data(0)[0],
        "image_0001.cbf",
    )

    assert imageset.get_format_class()

    imageset2 = ExperimentListFactory.from_filenames(["image_0001.cbf"])[0].imageset

    d_u_o = pytest.importorskip("dials.util.options")
    tolerance = d_u_o.tolerance_phil_scope.extract().tolerance

    print(sequence_diff(imageset, imageset2, tolerance=tolerance))

    assert BeamComparison()(imageset.get_beam(), imageset2.get_beam())
    assert DetectorComparison(origin_tolerance=tolerance.detector.origin)(
        imageset.get_detector(), imageset2.get_detector()
    )
    assert GoniometerComparison()(imageset.get_goniometer(), imageset2.get_goniometer())
    s1 = imageset.get_sequence()
    s2 = imageset.get_sequence()
    assert s1.get_exposure_times() == s2.get_exposure_times()
    assert s1.get_oscillation() == s2.get_oscillation()
    assert s1.get_image_range() == s2.get_image_range()
    assert imageset.get_raw_data(0) == imageset2.get_raw_data(0)
