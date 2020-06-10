from __future__ import absolute_import, division, print_function

import errno
import json
import os
from collections import namedtuple
from pprint import pprint

import mock
import pytest
import six.moves.cPickle as pickle

import dxtbx.datablock as datablock
from dxtbx.datablock import DataBlockFactory
from dxtbx.format.Format import Format
from dxtbx.imageset import ImageSequence
from dxtbx.model import Beam, Detector, Goniometer, Scan


@pytest.fixture(scope="session")
def centroid_test_data(dials_regression):
    return os.path.join(dials_regression, "centroid_test_data")


@pytest.fixture
def single_sequence_filenames(centroid_test_data):
    filenames = [
        os.path.join(centroid_test_data, "centroid_000{}.cbf".format(i))
        for i in range(1, 10)
    ]
    return filenames


@pytest.fixture
def multiple_sequence_filenames(centroid_test_data):
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
def multiple_block_filenames(single_sequence_filenames, all_image_examples):
    return single_sequence_filenames + all_image_examples


def encode_json_then_decode(obj, check_format=True):
    string = json.dumps([db.to_dict() for db in obj], ensure_ascii=True)
    return DataBlockFactory.from_json(string, check_format=check_format)


def pickle_then_unpickle(obj):
    return pickle.loads(pickle.dumps(obj))


def test_create_single_sequence(single_sequence_filenames):
    blocks = DataBlockFactory.from_filenames(single_sequence_filenames, verbose=True)
    assert len(blocks) == 1
    assert blocks[0].num_images() == 9
    assert blocks[0].format_class()
    imageset = blocks[0].extract_imagesets()
    assert len(imageset) == 1
    assert len(imageset[0]) == 9
    sequences = blocks[0].extract_sequences()
    assert len(sequences) == 1
    assert len(sequences[0]) == 9


def test_create_multiple_sequences(multiple_sequence_filenames):
    blocks = DataBlockFactory.from_filenames(multiple_sequence_filenames)
    assert len(blocks) == 1
    assert blocks[0].num_images() == 6
    assert blocks[0].format_class()
    imageset = blocks[0].extract_imagesets()
    assert len(imageset) == 2
    sequences = blocks[0].extract_sequences()
    assert len(sequences) == 2
    assert len(sequences[0]) == 3
    assert len(sequences[1]) == 3


def test_create_multiple_blocks(multiple_block_filenames):
    pprint(multiple_block_filenames)
    blocks = DataBlockFactory.from_filenames(multiple_block_filenames, verbose=True)
    assert blocks

    # Block 1
    assert blocks[0].num_images() == 12
    imageset = blocks[0].extract_imagesets()
    assert len(imageset) == 4
    assert len(imageset[0]) == 9
    sequences = blocks[0].extract_sequences()
    assert len(sequences) == 4
    assert len(sequences[0]) == 9

    pprint([b.num_images() for b in blocks])

    assert len(blocks) == 8


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
            if isinstance(im1, ImageSequence):
                assert isinstance(im2, ImageSequence)
                assert im1.get_beam() == im2.get_beam()
                assert im1.get_detector() == im2.get_detector()
                assert im1.get_goniometer() == im2.get_goniometer()
                assert im1.get_scan() == im2.get_scan()
            else:
                assert not isinstance(im2, ImageSequence)
                for i in range(len(im1)):
                    assert im1.get_beam(i) == im2.get_beam(i)
                    assert im1.get_detector(i) == im2.get_detector(i)


def test_from_null_sequence():
    filenames = ["template_%2d.cbf" % (i + 1) for i in range(0, 10)]
    sequence = Format.get_imageset(
        filenames,
        beam=Beam((0, 0, 1)),
        detector=Detector(),
        goniometer=Goniometer((1, 0, 0)),
        scan=Scan((1, 10), (0, 0.1)),
    )

    # Create the datablock
    datablock = DataBlockFactory.from_imageset(sequence)
    assert len(datablock) == 1
    datablock = datablock[0]
    assert datablock.format_class()

    sequences = datablock.extract_sequences()
    assert len(sequences) == 1
    assert sequences[0].get_beam() == sequence.get_beam()
    assert sequences[0].get_detector() == sequence.get_detector()
    assert sequences[0].get_goniometer() == sequence.get_goniometer()
    assert sequences[0].get_scan() == sequence.get_scan()


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


def test_path_iterator(monkeypatch):
    """Test the pathname iterator that avoids excessive file calls"""

    @classmethod
    def _fake_open_file(cls, name):
        """Mock replacement for Format's open_file"""
        if name in ["a", "b", "dir/c", "dir/d", "e"]:
            return mock.Mock()
        elif name.startswith("dir"):
            err = IOError()
            err.errno = errno.EISDIR
            # raise IOError(errno=errno.EISDIR)
            raise err
        assert False

    # Path the lookup of files
    listdir = mock.Mock(return_value=["c", "dir2", "d"])
    monkeypatch.setattr(os, "listdir", listdir)
    # Replace Format.open_file with a tame version
    monkeypatch.setattr(Format, "open_file", _fake_open_file)

    it = datablock._openingpathiterator(["a", "b", "dir", "e"])
    assert list(it) == ["a", "b", "dir/c", "dir/d", "e"]
    listdir.assert_called_once_with("dir")

    # Test that the list is sorted
    it = datablock._openingpathiterator(["e", "a", "b", "dir"])
    assert list(it) == ["a", "b", "dir/c", "dir/d", "e"]


def test_extract_metadata_record():
    """Make sure we can read a metadataobject from a format instance"""
    fmt = mock.MagicMock()
    fmt.get_image_file.return_value = "filename_000.cbf"
    fmt.get_scan.return_value = None
    record = datablock.ImageMetadataRecord.from_format(fmt)
    assert record.beam is fmt.get_beam()
    assert record.detector is fmt.get_detector()
    assert record.goniometer is fmt.get_goniometer()
    assert record.scan is None
    assert record.index is None


def _equal_but_not_same(thing):
    object_1 = tuple([thing])
    object_2 = tuple([thing])
    assert object_1 == object_2
    assert object_1 is not object_2
    return object_1, object_2


def test_merge_metadata_record():
    """Test that merging metadata records works correctly"""
    beam_a, beam_b = _equal_but_not_same("beam")
    detector_a, detector_b = _equal_but_not_same("detector")
    gonio_a, gonio_b = _equal_but_not_same("goniometer")

    a = datablock.ImageMetadataRecord(
        beam=beam_a, detector=detector_a, goniometer=gonio_a
    )
    b = datablock.ImageMetadataRecord(
        beam=beam_b, detector=detector_b, goniometer=gonio_b
    )
    pre_hash = hash(a)
    assert a.beam is not b.beam
    assert a.detector is not b.detector
    assert a.goniometer is not b.goniometer
    # This should do something
    assert b.merge_metadata_from(a)
    assert hash(a) == pre_hash, "a changed after merge"
    # Make sure metadata was merged
    assert a.beam is b.beam
    assert a.detector is b.detector
    assert a.goniometer is b.goniometer
    # This should NOT do something
    assert not a.merge_metadata_from(a)
    assert hash(a) == pre_hash


def test_merge_all_metadata():
    """Test that merging metadata over a whole list of records works"""
    beam_a, beam_b = _equal_but_not_same("beam")
    gonio_a, gonio_b = _equal_but_not_same("goniometer")

    a = datablock.ImageMetadataRecord(
        beam=beam_a, detector=object(), goniometer=gonio_a
    )
    b = datablock.ImageMetadataRecord(
        beam=beam_b, detector=object(), goniometer=gonio_b
    )
    records = [a, b]
    datablock._merge_model_metadata(records)
    assert a.beam is b.beam
    assert a.goniometer is b.goniometer
    assert a.detector is not b.detector


def test_merge_scan():
    """Test merging logic of scans"""
    scanA = mock.Mock(spec=Scan)
    scanB = mock.Mock(spec=Scan)
    recordA = mock.Mock(scan=scanA)
    recordB = mock.Mock(
        scan=scanB,
        beam=recordA.beam,
        detector=recordA.detector,
        goniometer=recordA.goniometer,
    )
    result = datablock._merge_scans([recordA, recordB])
    assert result == [recordA]
    scanA.append.assert_called_once_with(scanB)

    # Change some metadata in recordB so it doesn't match
    scanA.reset_mock()
    recordB.beam = mock.Mock()
    assert datablock._merge_scans([recordA, recordB]) == [recordA, recordB]


def test_groupby_template_none():
    Fo = namedtuple("Fo", ["template"])
    objs = [Fo(1), Fo(2), Fo(2), Fo(None), Fo(None), Fo("something")]
    result = list(datablock._groupby_template_is_none(objs))
    assert result == [
        [Fo(1)],
        [Fo(2)],
        [Fo(2)],
        [Fo(None), Fo(None)],
        [Fo("something")],
    ]


def test_iterate_with_previous():
    sample = list(range(5))
    assert list(datablock._iterate_with_previous(sample)) == [
        (None, 0),
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
    ]
