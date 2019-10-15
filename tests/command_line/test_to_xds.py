from __future__ import absolute_import, division, print_function

import procrunner
import pytest
from dxtbx.imageset import ImageSetFactory
from dxtbx.serialize import dump


@pytest.fixture(scope="session")
def expected_output(dials_data):
    return """\
DETECTOR=PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD=495976
SENSOR_THICKNESS= 0.320
DIRECTION_OF_DETECTOR_X-AXIS= 1.00000 0.00000 0.00000
DIRECTION_OF_DETECTOR_Y-AXIS= 0.00000 1.00000 0.00000
NX=2463 NY=2527 QX=0.1720 QY=0.1720
DETECTOR_DISTANCE= 190.180
ORGX= 1235.84 ORGY= 1279.58
ROTATION_AXIS= 1.00000 0.00000 0.00000
STARTING_ANGLE= 0.000
OSCILLATION_RANGE= 0.200
X-RAY_WAVELENGTH= 0.97950
INCIDENT_BEAM_DIRECTION= -0.000 -0.000 1.021
FRACTION_OF_POLARIZATION= 0.999
POLARIZATION_PLANE_NORMAL= 0.000 1.000 0.000
NAME_TEMPLATE_OF_DATA_FRAMES= %s
TRUSTED_REGION= 0.0 1.41
UNTRUSTED_RECTANGLE= 487 495 0 2528
UNTRUSTED_RECTANGLE= 981 989 0 2528
UNTRUSTED_RECTANGLE= 1475 1483 0 2528
UNTRUSTED_RECTANGLE= 1969 1977 0 2528
UNTRUSTED_RECTANGLE= 0 2464 195 213
UNTRUSTED_RECTANGLE= 0 2464 407 425
UNTRUSTED_RECTANGLE= 0 2464 619 637
UNTRUSTED_RECTANGLE= 0 2464 831 849
UNTRUSTED_RECTANGLE= 0 2464 1043 1061
UNTRUSTED_RECTANGLE= 0 2464 1255 1273
UNTRUSTED_RECTANGLE= 0 2464 1467 1485
UNTRUSTED_RECTANGLE= 0 2464 1679 1697
UNTRUSTED_RECTANGLE= 0 2464 1891 1909
UNTRUSTED_RECTANGLE= 0 2464 2103 2121
UNTRUSTED_RECTANGLE= 0 2464 2315 2333
DATA_RANGE= 1 9
JOB=XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT\
""" % (
        dials_data("centroid_test_data").join("centroid_????.cbf").strpath
    )


def test_to_xds_from_images(dials_data, expected_output, tmpdir):
    file_names = dials_data("centroid_test_data").listdir("centroid_*.cbf")
    result = procrunner.run(["dxtbx.to_xds"] + file_names, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    # allow extra lines to have been added (these may be comments)
    for record in expected_output.split("\n"):
        assert record.strip().encode("latin-1") in result.stdout, record


def test_to_xds_from_json(dials_data, expected_output, tmpdir):
    file_names = dials_data("centroid_test_data").listdir("centroid_*.cbf")

    # now test reading from a json file
    sequence = ImageSetFactory.new([f.strpath for f in file_names])[0]
    with tmpdir.join("sequence.json").open("wb") as fh:
        dump.imageset(sequence, fh)
    result = procrunner.run(["dxtbx.to_xds", "sequence.json"], working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    # allow extra lines to have been added (these may be comments)
    for record in expected_output.split("\n"):
        assert record.strip().encode("latin-1") in result.stdout, record
