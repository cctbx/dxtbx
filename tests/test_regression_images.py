"""
Image reading tests against the dials_regression suite
"""

import bz2
import gzip
import os
import shutil
from pathlib import Path

import pytest

import scitbx.matrix
from rstbx.slip_viewer.slip_viewer_image_factory import SlipViewerImageFactory
from scitbx.array_family import flex

import dxtbx.format.Registry
from dxtbx.util import get_url_scheme

try:
    import h5py
except ModuleNotFoundError:
    h5py = None

try:
    import xfel
except ImportError:
    # Although ModuleNotFoundError is correct, cctbx may raise an ImportError
    # when the xfel extension module is not found.
    xfel = None

try:
    import cbflib_adaptbx
except ModuleNotFoundError:
    # Detectorbase for CBF won't work
    cbflib_adaptbx = None

_files = (
    "ALS_1231/q315r_lyso_1_001.img",
    "ALS_422/lyso_041013a_1_001.img",
    "ALS_501/als501_q4_1_001.img",
    "ALS_733/200mMNaCl5pcGlyc_400.edf",
    "ALS_821/q210_lyso_1_101.img",
    "ALS_831/q315r_lyso_001.img",
    "APS_14BMC/q315_1_001.img",
    "APS_17ID/q210_1_001.img",
    "APS_19ID/q315_unbinned_a.0001.img",
    "APS_19ID/t1.0001.img.bz2",
    "APS_22ID/mar300.0001",
    "APS_23IDD/mar300_1_E1.0001",
    "APS_24IDC/pilatus_1_0001.cbf",
    "APS_24IDC/q315_1_001.img",
    "APS_24IDE_test/thaum-12_1_0001.cbf",
    "APS_24IDE_test/thaum-12_1_0002.cbf",
    "APS_24IDE_test/thaum-12_1_0003.cbf",
    "APS_24IDE_test/thaum-12_1_0004.cbf",
    "APS_24IDE_test/thaum-12_1_0005.cbf",
    "APS_24IDE_test/thaum-12_1_0006.cbf",
    "APS_24IDE_test/thaum-12_1_0007.cbf",
    "APS_24IDE_test/thaum-12_1_0008.cbf",
    "APS_24IDE_test/thaum-12_1_0009.cbf",
    "APS_24IDE_test/thaum-12_1_0010.cbf",
    "Bruker_PHOTON_II/dan_01_0001.sfrm",
    "CLS1_08ID1/mar225_2_E0_0001.img",
    "DESY_BW7B/mar345_01_001.mar2300",
    "DESY_ID141/q210_2_001.img",
    "DLS_I02/X4_wide_M1S4_1_0001.cbf",
    "DLS_I04/grid_full_cbf_0005.cbf",
    "DLS_I19/I19_P300k_00001.cbf",
    "DLS_I23/I23_P12M_alpha_0001.cbf",
    "DLS_I23/germ_13KeV_0001.cbf",
    "DLS_I24_stills/still_0001.cbf",
    "DLS_eBIC/nanB_4_singleframe_1.ser.bz2",
    "DTREK_home_lab/s01f0001.osc",
    "ESRF_BM14/mar165_001.mccd",
    "ESRF_BM14/mar225_1_001.mccd",
    "ESRF_ID231/q315r_7_001.img",
    "ESRF_ID29/trypsin_1_0001.cbf",
    "LCLS_cspad_nexus/idx-20130301060858401.cbf",
    "LCLS_cspad_nexus/idx-20130301060858601.cbf",
    "LCLS_cspad_nexus/idx-20130301060858701.cbf",
    "LCLS_cspad_nexus/idx-20130301060858801.cbf",
    "LCLS_jungfrau/jungfrau_multipanel.cbf",
    "MLFSOM_simulation/fake_00001.img",
    "MacScience/reallysurprise_001.ipf",
    "RAXIS-HTC/test1_lysozyme_0111060001.osc",
    "RigakuA200/XV10001.img",
    "SACLA_MPCCD_Cheetah/run266702-0-subset.h5",
    "SLS_Eiger_16M_as_CBF/insu_with_bs_labelit_0001.cbf",
    "SLS_Eiger_16M_as_CBF/insu_with_bs_labelit_0901.cbf",
    "SLS_X06SA/mar225_2_001.img",
    "SLS_X06SA/pilatus6m_1_00001.cbf",
    "SPring8_ADSC_SN916/Xtal17-2phi_3_015.cbf",
    "SPring8_BL12B2_MX225HE/lys001_000001.img",
    "SPring8_BL12B2_MX225HE/lys001_000091.img",
    "SPring8_BL26B1_Raxis5/raxis5_000001.img",
    "SPring8_BL26B1_Raxis5/raxis5_000091.img",
    "SPring8_BL26B1_SaturnA200/A200_000001.img",
    "SPring8_BL26B1_SaturnA200/A200_000002.img",
    "SPring8_BL26B2_MX225/2sec_Al200um_000001.img",
    "SPring8_BL26B2_MX225/2sec_Al200um_000090.img",
    "SPring8_BL32XU/lys_00001.img",
    "SPring8_BL32XU/rayonix225_0001.img",
    "SPring8_BL32XU/rayonix225hs_0001.img",
    "SPring8_BL32XU_MX225HS/ds_000001.img",
    "SPring8_BL32XU_MX225HS/ds_000045.img",
    "SPring8_BL38B1_MX225HE/bl38b1_001.img",
    "SPring8_BL38B1_MX225HE/bl38b1_090.img",
    "SPring8_BL41XU_PILATUS3_6M/data1_000001.cbf",
    "SPring8_BL41XU_PILATUS3_6M/data1_000901.cbf",
    "SPring8_BL44XU_MX300HE/bl44xu_lys_000001.img",
    "SPring8_BL44XU_MX300HE/bl44xu_lys_000002.img",
    "SRS_101/mar225_001.img",
    "SRS_142/q4_1_001.img",
    "SSRL_bl111/mar325_1_001.mccd",
    "SSRL_bl91/q315_1_001.img",
    "Texas_A_and_M_University/lyziph6p5_01_0001.sfrm",
    "XDS/INTEGRATE.HKL",
    "XDS/XDS_ASCII.HKL",
    "XDS/XPARM.XDS",
    "dials-190/whatev1_01_00001.cbf",
    "dials-190/whatev1_01_00002.cbf",
    "dials-190/whatev1_02_00001.cbf",
    "dials-190/whatev1_02_00002.cbf",
    "dials-190/whatev1_03_00001.cbf",
    "dials-190/whatev1_03_00002.cbf",
    "johns_hopkins_raxisII/lys_001.osc",
    "saturn/lyso_00001.img",
    "xia2/merge2cbf_averaged_0001.cbf",
)

_files_with_detectorbase = sorted(
    set(_files)
    - {
        "Bruker_PHOTON_II/dan_01_0001.sfrm",
        "DLS_eBIC/nanB_4_singleframe_1.ser.bz2",
        "LCLS_cspad_nexus/idx-20130301060858401.cbf",
        "LCLS_cspad_nexus/idx-20130301060858601.cbf",
        "LCLS_cspad_nexus/idx-20130301060858701.cbf",
        "LCLS_cspad_nexus/idx-20130301060858801.cbf",
        "LCLS_jungfrau/jungfrau_multipanel.cbf",
        "SACLA_MPCCD_Cheetah/run266702-0-subset.h5",
        "XDS/INTEGRATE.HKL",
        "XDS/XDS_ASCII.HKL",
        "XDS/XPARM.XDS",
    }
)


def test_berkeley_special_h5():
    # Handle the special berkeley-only h5 file
    special_h5 = "/net/viper/raid1/dectris/eiger16MNov2015/2015_11_10/insu6_1_master.h5"
    if not os.path.isfile(special_h5):
        pytest.skip("LBL-only file not present")

    # Run the tests, but without dials_regression fixture, circumventing the path resolution step
    test_read_image(special_h5, dials_regression=None)
    test_format_class_API_assumptions(special_h5, dials_regression=None)


@pytest.mark.regression
@pytest.mark.parametrize("test_image", _files_with_detectorbase)
def test_detectorbase(test_image, dials_regression):
    if not xfel and test_image.startswith("LCLS"):
        pytest.skip("could not import 'xfel'")

    if not h5py and test_image.endswith((".h5", ".nxs")):
        pytest.skip("could not import 'h5py'")

    if test_image.endswith(".cbf") and not cbflib_adaptbx:
        pytest.skip("No cbflib_adaptbx: CBF Detectorbase is not available")

    test_image = str(
        Path(dials_regression).joinpath("image_examples", *test_image.split("/"))
    )
    format_instance = dxtbx.format.Registry.get_format_class_for_file(test_image)
    print("Reading", test_image)
    print("Format:", format_instance)
    assert format_instance, "no matching format class found"
    instance = format_instance(test_image)

    instance.get_detectorbase()

    imgfactory = SlipViewerImageFactory(test_image)
    imgfactory.read()
    print("  Detectorbase:", instance.detectorbase.__class__.__name__)

    try:
        print(imgfactory.rawdata.focus())
    except AttributeError:
        # Not all instances have this attribute
        print("  multireadout")

    I_raw_data = imgfactory.get_raw_data()
    if not isinstance(I_raw_data, tuple):
        I_raw_data = (I_raw_data,)

    try:
        R_raw_data = instance.get_raw_data()
    except TypeError:
        R_raw_data = instance.get_raw_data(0)

    if not isinstance(R_raw_data, tuple):
        R_raw_data = (R_raw_data,)

    # NOTE dxtbx and image factory arrays are compared here for identical values.
    for Ip, Rp in zip(I_raw_data, R_raw_data):
        assert (Ip == Rp).all_eq(True)


@pytest.mark.regression
@pytest.mark.parametrize("test_image", _files)
def test_read_image(test_image, dials_regression):
    """Test that all the regression images can be read"""
    if not xfel and test_image.startswith("LCLS"):
        pytest.skip("could not import 'xfel'")

    if not h5py and test_image.endswith((".h5", ".nxs")):
        pytest.skip("could not import 'h5py'")

    if dials_regression:
        test_image = os.path.join(
            dials_regression, "image_examples", *test_image.split("/")
        )

    format_instance = dxtbx.format.Registry.get_format_class_for_file(test_image)
    print("Reading", test_image)
    print("Format:", format_instance)
    assert format_instance, "no matching format class found"
    instance = format_instance(test_image)

    # Test metadata reading
    instance.get_goniometer()
    instance.get_beam()
    instance.get_scan()
    detector = instance.get_detector()
    # From old test_dxtbx; get the image size
    if detector is not None:
        detector[0].get_image_size()

    for panel in detector:
        d_mat = scitbx.matrix.sqr(panel.get_d_matrix())
        if d_mat.determinant() <= 0:
            print("  d matrix with non-positive determinant")

    # Test reading of the raw data
    # XDS, HKL we expect to fail for this  - so skip this part for those
    if not test_image.endswith(("XDS", "HKL")):
        try:
            R_raw_data = instance.get_raw_data()
        except TypeError:
            R_raw_data = instance.get_raw_data(0)

        if not isinstance(R_raw_data, tuple):
            R_raw_data = (R_raw_data,)

        print("%-40s" % format_instance.__name__, R_raw_data[0].focus())

        # Specific test for cctbx/dxtbx#163. This test will fail if char is unsigned.
        if "APS_24IDC" in test_image and "pilatus_1_0001.cbf" in test_image:
            d = R_raw_data[0]
            assert (
                flex.sum(d.as_1d().select(d.as_1d() >= 0)) == 20108255
            )  # matches total counts from dxtbx.print_header


@pytest.mark.parametrize("test_image", _files)
def test_format_class_API_assumptions(test_image, dials_regression):
    """For a given image file, walk the whole DAG of Format objects, and
    verify the following basic assumptions for format classes:
    * Any file must only be understood by a single leaf node format class.
    * No .understand() call on any top level format class or a child class
      of another understanding format is allowed to throw an exception.
    """
    if not h5py and test_image.endswith((".h5", ".nxs")):
        pytest.skip("could not import 'h5py'")

    if dials_regression:
        test_image = os.path.join(
            dials_regression, "image_examples", *test_image.split("/")
        )

    dag = dxtbx.format.Registry.get_format_class_dag()

    def recurse(parentformat, filename, level=0):
        known_format_class = None
        multiple_formats = False
        for subformat in dag.get(parentformat, []):
            format_class = dxtbx.format.Registry.get_format_class_for(subformat)
            if not get_url_scheme(filename) in format_class.schemes:
                print("Not matching ", filename, "to", format_class)
                continue
            understood = format_class.understand(filename)
            print("%s%s: %s" % ("  " * level, subformat, understood))
            if understood:
                recursive_format_class, subtree_multiple = recurse(
                    subformat, filename, level + 1
                )
                understood_format_class = recursive_format_class or subformat
                if known_format_class and known_format_class != understood_format_class:
                    print(
                        f"File can be understood as {known_format_class} and {understood_format_class}"
                    )
                    multiple_formats = True
                known_format_class = understood_format_class
                multiple_formats |= subtree_multiple
        return known_format_class, multiple_formats

    understood_format, multiple_formats = recurse("Format", test_image)

    assert not multiple_formats, "image file understood by multiple Format objects"
    # It's a failure if nothing could understand this file
    assert understood_format, "No formatter could be found"
    print("File understood as", understood_format)


@pytest.mark.regression
@pytest.mark.parametrize(
    "image_to_test",
    # Test miniCBF and a flavour of FullCBF
    ["DLS_I02/X4_wide_M1S4_1_0001.cbf", "DLS_I04/grid_full_cbf_0005.cbf"],
)
@pytest.mark.parametrize(
    "compression,extension",
    [(gzip.GzipFile, "gz"), (bz2.BZ2File, "bz2")],
    ids=["gzip", "bz2"],
)
def test_compressed_images(
    dials_regression, compression, extension, image_to_test, tmp_path
):
    path = Path(dials_regression) / "image_examples" / image_to_test
    if not path.exists():
        pytest.skip(str(path) + " not present in dials_regression")

    # Compress this file and write to disk
    test_image_for_reading = str(tmp_path / (path.name + "." + extension))
    with path.open("rb") as f_in:
        with compression(test_image_for_reading, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Find and instantiate the format class
    format_class = dxtbx.format.Registry.get_format_class_for_file(
        test_image_for_reading
    )
    assert format_class is not None, "no matching format class found"
    instance = format_class(test_image_for_reading)

    # Test metadata reading
    assert instance.get_goniometer()
    assert instance.get_beam()
    assert instance.get_scan()
    assert instance.get_detector()
    # Test data reading
    assert instance.get_raw_data()
