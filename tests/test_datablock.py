from __future__ import absolute_import, division, print_function

import json
import os
from builtins import range
from pprint import pprint

import pytest
import six.moves.cPickle as pickle

from dxtbx.datablock import DataBlockFactory
from dxtbx.format.Format import Format
from dxtbx.imageset import ImageSweep
from dxtbx.model import Beam, Detector, Goniometer, Scan


@pytest.fixture(scope="session")
def centroid_test_data(dials_regression):
    return os.path.join(dials_regression, "centroid_test_data")


@pytest.fixture
def single_sweep_filenames(centroid_test_data):
    filenames = [
        os.path.join(centroid_test_data, "centroid_000{}.cbf".format(i))
        for i in range(1, 10)
    ]
    return filenames


@pytest.fixture
def multiple_sweep_filenames(centroid_test_data):
    filenames = [
        os.path.join(centroid_test_data, "centroid_000{}.cbf".format(i))
        for i in [1, 2, 3, 7, 8, 9]
    ]
    return filenames


@pytest.fixture
def all_image_examples(dials_regression):
    filenames = (
        ("ALS_1231", "q315r_lyso_1_001.img"),
        ("ALS_501", "als501_q4_1_001.img"),
        ("ALS_821", "q210_lyso_1_101.img"),
        ("ALS_831", "q315r_lyso_001.img"),
        ("APS_14BMC", "q315_1_001.img"),
        ("APS_17ID", "q210_1_001.img"),
        ("APS_19ID", "q315_unbinned_a.0001.img"),
        ("APS_22ID", "mar300.0001"),
        ("APS_23IDD", "mar300_1_E1.0001"),
        ("APS_24IDC", "pilatus_1_0001.cbf"),
        ("APS_24IDC", "q315_1_001.img"),
        ("CLS1_08ID1", "mar225_2_E0_0001.img"),
        ("DESY_ID141", "q210_2_001.img"),
        ("ESRF_BM14", "mar165_001.mccd"),
        ("ESRF_BM14", "mar225_1_001.mccd"),
        ("ESRF_ID231", "q315r_7_001.img"),
        ("RAXIS-HTC", "test1_lysozyme_0111060001.osc"),
        ("SLS_X06SA", "mar225_2_001.img"),
        ("SLS_X06SA", "pilatus6m_1_00001.cbf"),
        ("SRS_101", "mar225_001.img"),
        ("SRS_142", "q4_1_001.img"),
        ("SSRL_bl111", "mar325_1_001.mccd"),
        ("xia2", "merge2cbf_averaged_0001.cbf"),
    )
    return [os.path.join(dials_regression, "image_examples", *f) for f in filenames]


@pytest.fixture
def multiple_block_filenames(single_sweep_filenames, all_image_examples):
    return single_sweep_filenames + all_image_examples


def encode_json_then_decode(obj, check_format=True):
    string = json.dumps([db.to_dict() for db in obj], ensure_ascii=True)
    return DataBlockFactory.from_json(string, check_format=check_format)


def pickle_then_unpickle(obj):
    return pickle.loads(pickle.dumps(obj))


def test_create_single_sweep(single_sweep_filenames):
    blocks = DataBlockFactory.from_filenames(single_sweep_filenames, verbose=True)
    assert len(blocks) == 1
    assert blocks[0].num_images() == 9
    assert blocks[0].format_class()
    imageset = blocks[0].extract_imagesets()
    assert len(imageset) == 1
    assert len(imageset[0]) == 9
    sweeps = blocks[0].extract_sweeps()
    assert len(sweeps) == 1
    assert len(sweeps[0]) == 9


def test_create_multiple_sweeps(multiple_sweep_filenames):
    blocks = DataBlockFactory.from_filenames(multiple_sweep_filenames)
    assert len(blocks) == 1
    assert blocks[0].num_images() == 6
    assert blocks[0].format_class()
    imageset = blocks[0].extract_imagesets()
    assert len(imageset) == 2
    sweeps = blocks[0].extract_sweeps()
    assert len(sweeps) == 2
    assert len(sweeps[0]) == 3
    assert len(sweeps[1]) == 3


def test_create_multiple_blocks(multiple_block_filenames):
    pprint(multiple_block_filenames)
    blocks = DataBlockFactory.from_filenames(multiple_block_filenames, verbose=True)
    assert blocks

    # Block 1
    assert blocks[0].num_images() == 9
    imageset = blocks[0].extract_imagesets()
    assert len(imageset) == 1
    assert len(imageset[0]) == 9
    sweeps = blocks[0].extract_sweeps()
    assert len(sweeps) == 1
    assert len(sweeps[0]) == 9

    pprint([b.num_images() for b in blocks])

    assert len(blocks) == 22


def test_pickling(multiple_block_filenames):
    blocks1 = DataBlockFactory.from_filenames(multiple_block_filenames, verbose=True)
    blocks2 = pickle_then_unpickle(blocks1)
    assert len(blocks2) == len(blocks1)
    for b1, b2 in zip(blocks1, blocks2):
        assert b1.format_class() == b2.format_class()
        assert b1 == b2
    assert blocks1 == blocks2


def test_json(multiple_block_filenames):
    blocks1 = DataBlockFactory.from_filenames(multiple_block_filenames, verbose=True)
    blocks2 = encode_json_then_decode(blocks1)
    assert len(blocks2) == len(blocks1)
    for b1, b2 in zip(blocks1, blocks2):
        assert b1.format_class() == b2.format_class()
        assert b1 == b2
    assert blocks1 == blocks2


def test_json2(multiple_block_filenames):
    blocks1 = DataBlockFactory.from_filenames(multiple_block_filenames, verbose=True)
    blocks2 = encode_json_then_decode(blocks1, check_format=False)
    assert len(blocks2) == len(blocks1)
    for b1, b2 in zip(blocks1, blocks2):
        for im1, im2 in zip(b1.extract_imagesets(), b2.extract_imagesets()):
            assert len(im1) == len(im2)
            if isinstance(im1, ImageSweep):
                assert isinstance(im2, ImageSweep)
                assert im1.get_beam() == im2.get_beam()
                assert im1.get_detector() == im2.get_detector()
                assert im1.get_goniometer() == im2.get_goniometer()
                assert im1.get_scan() == im2.get_scan()
            else:
                assert not isinstance(im2, ImageSweep)
                for i in range(len(im1)):
                    assert im1.get_beam(i) == im2.get_beam(i)
                    assert im1.get_detector(i) == im2.get_detector(i)


def test_from_null_sweep():
    filenames = ["template_%2d.cbf" % (i + 1) for i in range(0, 10)]
    sweep = Format.get_imageset(
        filenames,
        beam=Beam((0, 0, 1)),
        detector=Detector(),
        goniometer=Goniometer((1, 0, 0)),
        scan=Scan((1, 10), (0, 0.1)),
    )

    # Create the datablock
    datablock = DataBlockFactory.from_imageset(sweep)
    assert len(datablock) == 1
    datablock = datablock[0]
    assert datablock.format_class()

    sweeps = datablock.extract_sweeps()
    assert len(sweeps) == 1
    assert sweeps[0].get_beam() == sweep.get_beam()
    assert sweeps[0].get_detector() == sweep.get_detector()
    assert sweeps[0].get_goniometer() == sweep.get_goniometer()
    assert sweeps[0].get_scan() == sweep.get_scan()


def test_with_bad_external_lookup(centroid_test_data):
    filename = os.path.join(centroid_test_data, "datablock_with_bad_lookup.json")
    blocks = DataBlockFactory.from_json_file(filename, check_format=False)
    assert len(blocks) == 1
    imageset = blocks[0].extract_imagesets()[0]
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.empty()
    assert imageset.external_lookup.gain.data.empty()
    assert imageset.external_lookup.pedestal.data.empty()

    blocks = encode_json_then_decode(blocks, check_format=False)
    assert len(blocks) == 1
    imageset = blocks[0].extract_imagesets()[0]
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.empty()
    assert imageset.external_lookup.gain.data.empty()
    assert imageset.external_lookup.pedestal.data.empty()


def test_with_external_lookup(centroid_test_data):
    filename = os.path.join(centroid_test_data, "datablock_with_lookup.json")
    blocks = DataBlockFactory.from_json_file(filename)
    assert len(blocks) == 1
    imageset = blocks[0].extract_imagesets()[0]
    assert not imageset.external_lookup.mask.data.empty()
    assert not imageset.external_lookup.gain.data.empty()
    assert not imageset.external_lookup.pedestal.data.empty()
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.tile(0).data().all_eq(True)
    assert imageset.external_lookup.gain.data.tile(0).data().all_eq(1)
    assert imageset.external_lookup.pedestal.data.tile(0).data().all_eq(0)

    blocks = encode_json_then_decode(blocks)
    assert len(blocks) == 1
    imageset = blocks[0].extract_imagesets()[0]
    assert not imageset.external_lookup.mask.data.empty()
    assert not imageset.external_lookup.gain.data.empty()
    assert not imageset.external_lookup.pedestal.data.empty()
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.tile(0).data().all_eq(True)
    assert imageset.external_lookup.gain.data.tile(0).data().all_eq(1)
    assert imageset.external_lookup.pedestal.data.tile(0).data().all_eq(0)
