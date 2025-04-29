"""
Image reading tests against image_examples data
"""

from __future__ import annotations

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
    import cbflib_adaptbx

    # Check this isn't a namespace package
    if hasattr(cbflib_adaptbx, "__path__") and cbflib_adaptbx.__file__ is None:
        cbflib_adaptbx = None
except ModuleNotFoundError:
    # Detectorbase for CBF won't work
    cbflib_adaptbx = None

from . import imagelist

_files = imagelist.image_examples


@pytest.mark.regression
@pytest.mark.parametrize("test_image", _files)
def test_detectorbase(test_image):
    if not h5py and test_image.endswith((".h5", ".nxs")):
        pytest.skip("could not import 'h5py'")

    if test_image.endswith(".cbf") and not cbflib_adaptbx:
        pytest.skip("No cbflib_adaptbx: CBF Detectorbase is not available")

    if test_image.endswith("bz2") or test_image.endswith("gz"):
        # Skip compressed files
        pytest.skip("Skip compressed files as detectorbase cannot be instantiated")

    format_instance = dxtbx.format.Registry.get_format_class_for_file(test_image)
    print("Reading", test_image)
    print("Format:", format_instance)
    assert format_instance, "no matching format class found"
    instance = format_instance(test_image)

    try:
        detectorbase = instance.get_detectorbase()
    except NotImplementedError:
        pytest.skip("Skipping image with no detectorbase")

    if not detectorbase:
        pytest.skip("Skipping image with no detectorbase")

    imgfactory = SlipViewerImageFactory(test_image)
    try:
        imgfactory.read()
    except AttributeError:
        pytest.skip("Skip case where SlipViewerImageFactory cannot read")

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
def test_read_image(test_image):
    """Test that all the regression images can be read"""

    if not h5py and test_image.endswith((".h5", ".nxs")):
        pytest.skip("could not import 'h5py'")

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
def test_format_class_API_assumptions(test_image):
    """For a given image file, walk the whole DAG of Format objects, and
    verify the following basic assumptions for format classes:
    * Any file must only be understood by a single leaf node format class.
    * No .understand() call on any top level format class or a child class
      of another understanding format is allowed to throw an exception.
    """
    if not h5py and test_image.endswith((".h5", ".nxs")):
        pytest.skip("could not import 'h5py'")

    dag = dxtbx.format.Registry.get_format_class_dag()

    def recurse(parentformat, filename, level=0):
        known_format_class = None
        multiple_formats = False
        for subformat in dag.get(parentformat, []):
            format_class = dxtbx.format.Registry.get_format_class_for(subformat)
            if get_url_scheme(filename) not in format_class.schemes:
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
