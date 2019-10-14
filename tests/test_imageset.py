from __future__ import absolute_import, division, print_function

import os
from builtins import range

from scitbx.array_family import flex

import dxtbx.format.image
import dxtbx.format.Registry
import dxtbx.tests.imagelist
import pytest
import six.moves.cPickle as pickle
from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus as FormatClass
from dxtbx.imageset import ExternalLookup, ImageSequence, ImageSetData, ImageSetFactory
from dxtbx.model import Beam, Detector, Panel


@pytest.mark.parametrize(
    "image",
    dxtbx.tests.imagelist.smv_images
    + dxtbx.tests.imagelist.tiff_images
    + dxtbx.tests.imagelist.cbf_multitile_images
    + dxtbx.tests.imagelist.cbf_images,
    ids=(
        dxtbx.tests.imagelist.smv_image_ids
        + dxtbx.tests.imagelist.tiff_image_ids
        + dxtbx.tests.imagelist.cbf_multitile_image_ids
        + dxtbx.tests.imagelist.cbf_image_ids
    ),
)
def test_format(dials_regression, image):
    print(image)
    image = os.path.join(dials_regression, *(image.split("/")))
    format_class = dxtbx.format.Registry.get_format_class_for_file(image)
    reader = format_class.get_reader()([image])

    N = len(reader)

    for i in range(N):
        reader.read(i)

    assert format_class.get_imageset([image])


def test_image_tile():
    data = flex.int(flex.grid(10, 10))
    name = "TileName"

    tile = dxtbx.format.image.ImageTileInt(data, name)

    assert tile.data().all_eq(data)
    assert tile.name() == name
    assert tile.empty() is False


def test_image():
    data = flex.int(flex.grid(10, 10))
    name = "TileName0"
    tile0 = dxtbx.format.image.ImageTileInt(data, name)
    image = dxtbx.format.image.ImageInt(tile0)
    for i in range(1, 4):
        data = flex.int(flex.grid(10, 10))
        name = "TileName%d" % i
        tile = dxtbx.format.image.ImageTileInt(data, name)
        image.append(tile)

    assert image.n_tiles() == 4
    for i in range(image.n_tiles()):
        tile = image.tile(i)
        assert tile.name() == "TileName%d" % i


def test_image_buffer():
    data = flex.int(flex.grid(10, 10))
    name = "TileName0"
    tile0 = dxtbx.format.image.ImageTileInt(data, name)
    image = dxtbx.format.image.ImageInt(tile0)

    b = dxtbx.format.image.ImageBuffer(image)
    assert b.is_int() is True
    assert b.is_double() is False
    assert b.is_empty() is False


def test_external_lookup():
    mask = flex.bool(flex.grid(10, 10), True)
    gain = flex.double(flex.grid(10, 10), 1)
    pedestal = flex.double(flex.grid(10, 10), 2)

    lookup = ExternalLookup()
    lookup.mask.data = dxtbx.format.image.ImageBool(
        dxtbx.format.image.ImageTileBool(mask)
    )
    lookup.gain.data = dxtbx.format.image.ImageDouble(
        dxtbx.format.image.ImageTileDouble(gain)
    )
    lookup.pedestal.data = dxtbx.format.image.ImageDouble(
        dxtbx.format.image.ImageTileDouble(pedestal)
    )

    mask2 = lookup.mask.data.tile(0).data()
    gain2 = lookup.gain.data.tile(0).data()
    pedestal2 = lookup.pedestal.data.tile(0).data()

    assert mask2.all_eq(mask)
    assert gain2.all_eq(gain)
    assert pedestal2.all_eq(pedestal)


def test_imagesetdata(centroid_files):
    ReaderClass = FormatClass.get_reader()

    reader = ReaderClass(centroid_files)
    masker = FormatClass(centroid_files[0]).get_masker()

    handle = ImageSetData(reader, masker)

    assert handle.get_data(0).as_int().tile(0).data()

    assert handle.has_single_file_reader() is False

    path = handle.get_path(0)
    assert path == centroid_files[0]

    master_path = handle.get_master_path()
    assert master_path == ""

    identifier = handle.get_image_identifier(0)
    assert identifier == centroid_files[0]

    beam = FormatClass(centroid_files[0]).get_beam()
    detector = FormatClass(centroid_files[0]).get_detector()
    goniometer = FormatClass(centroid_files[0]).get_goniometer()
    scan = FormatClass(centroid_files[0]).get_scan()

    handle.set_beam(beam, 0)
    handle.set_detector(detector, 0)
    handle.set_goniometer(goniometer, 0)
    handle.set_scan(scan, 0)

    beam2 = handle.get_beam(0)
    detector2 = handle.get_detector(0)
    goniometer2 = handle.get_goniometer(0)
    scan2 = handle.get_scan(0)

    assert beam2 == beam
    assert detector2 == detector
    assert goniometer2 == goniometer
    assert scan2 == scan

    mask = flex.bool(flex.grid(10, 10), True)
    gain = flex.double(flex.grid(10, 10), 1)
    pedestal = flex.double(flex.grid(10, 10), 2)

    handle.external_lookup.mask.data = dxtbx.format.image.ImageBool(
        dxtbx.format.image.ImageTileBool(mask)
    )
    handle.external_lookup.gain.data = dxtbx.format.image.ImageDouble(
        dxtbx.format.image.ImageTileDouble(gain)
    )
    handle.external_lookup.pedestal.data = dxtbx.format.image.ImageDouble(
        dxtbx.format.image.ImageTileDouble(pedestal)
    )

    mask2 = handle.external_lookup.mask.data.tile(0).data()
    gain2 = handle.external_lookup.gain.data.tile(0).data()
    pedestal2 = handle.external_lookup.pedestal.data.tile(0).data()

    assert mask2.all_eq(mask)
    assert gain2.all_eq(gain)
    assert pedestal2.all_eq(pedestal)


@pytest.fixture(scope="session")
def centroid_files(dials_data):
    return [
        dials_data("centroid_test_data").join("centroid_%04d.cbf" % i).strpath
        for i in range(1, 10)
    ]


@pytest.fixture
def centroid_files_and_imageset(centroid_files):
    # Create the format class
    format_class = dxtbx.format.Registry.get_format_class_for_file(centroid_files[0])

    # Create the reader
    imageset = format_class.get_imageset(centroid_files, as_imageset=True)

    return centroid_files, imageset


def assert_is_iterable(iterator):
    list(iterator)


def assert_can_get_detectorbase(obj, indices, outside_index):
    for i in indices:
        obj.get_detectorbase(i)

    with pytest.raises(RuntimeError):
        obj.get_detectorbase(outside_index)


class TestImageSet(object):
    def test_imageset(self, centroid_files_and_imageset):
        filenames, imageset = centroid_files_and_imageset

        # Run a load of tests
        self.tst_get_item(imageset)
        assert len(imageset) == len(filenames)
        assert_is_iterable(imageset)
        self.tst_paths(imageset, filenames)
        assert_can_get_detectorbase(imageset, range(len(filenames)), 9)
        self.tst_get_models(imageset, list(range(len(filenames))), 9)

    def tst_get_item(self, imageset):
        assert imageset[0]
        with pytest.raises(RuntimeError):
            _ = imageset[9]

        imageset2 = imageset[3:7]
        assert imageset2[0]
        with pytest.raises(RuntimeError):
            _ = imageset2[5]

        assert len(imageset2) == 4
        assert_can_get_detectorbase(imageset2, range(0, 4), 5)
        self.tst_get_models(imageset2, range(0, 4), 5)
        self.tst_paths(imageset2, imageset.paths()[3:7])
        assert_is_iterable(imageset2)

        imageset2 = imageset[3:5]
        assert imageset2[0]
        with pytest.raises(RuntimeError):
            _ = imageset2[2]

        assert len(imageset2) == 2
        assert_can_get_detectorbase(imageset2, range(0, 2), 2)
        self.tst_get_models(imageset2, range(0, 2), 2)
        self.tst_paths(imageset2, imageset.paths()[3:5])
        assert_is_iterable(imageset2)

    @staticmethod
    def tst_paths(imageset, filenames1):
        filenames2 = imageset.paths()
        for f1, f2 in zip(filenames1, filenames2):
            assert f1 == f2

    def tst_get_models(self, imageset, indices, outside_index):
        for i in indices:
            self.tst_get_models_index(imageset, i)

        with pytest.raises(RuntimeError):
            self.tst_get_models_index(imageset, outside_index)

    @staticmethod
    def tst_get_models_index(imageset, index=None):
        imageset.get_detector(index)
        imageset.get_beam(index)

    @staticmethod
    def tst_set_models(imageset):
        # Create some other models
        beam = Beam((1, 0, 0), 0.5)
        detector = Detector(
            Panel(
                "UNKNOWN",
                "Panel",
                (1, 0, 0),
                (0, 1, 0),
                (0, 0, 1),
                (0.1, 0.1),
                (1000, 1000),
                (0, 1),
            )
        )

        # Override sequence models
        imageset.set_beam(beam)
        imageset.set_detector(detector)

        # Ensure this doens't interfere with reading
        for i in imageset:
            pass

        # Get the models back and check they're ok
        beam2 = imageset.get_beam()
        detector2 = imageset.get_detector()
        assert beam2 == beam
        assert detector2 == detector

        # Get the models from an index back and check they're not the same
        beam2 = imageset.get_beam(0)
        detector2 = imageset.get_detector(0)
        assert beam2 != beam
        assert detector2 != detector


class TestImageSequence(object):
    def test(self, centroid_files):
        # Create the format class
        format_class = dxtbx.format.Registry.get_format_class_for_file(
            centroid_files[0]
        )

        # Create the sequence
        sequence = format_class.get_imageset(centroid_files)

        # Run a load of tests
        assert len(sequence) == len(centroid_files)
        assert sequence.get_array_range() == (0, 9)
        self.tst_get_item(sequence)
        assert_is_iterable(sequence)
        self.tst_paths(sequence, centroid_files)
        assert_can_get_detectorbase(sequence, range(len(centroid_files)), 9)
        self.tst_get_models(sequence, range(len(centroid_files)), 9)
        self.tst_set_models(sequence)

    def tst_get_item(self, sequence):
        _ = sequence[0]
        with pytest.raises(RuntimeError):
            _ = sequence[9]

        sequence2 = sequence[3:7]
        assert sequence2.get_array_range() == (3, 7)
        _ = sequence2[0]
        with pytest.raises(RuntimeError):
            _ = sequence2[5]

        assert len(sequence2) == 4
        assert_can_get_detectorbase(sequence2, range(0, 4), 5)
        self.tst_get_models(sequence2, range(0, 4), 5)
        self.tst_paths(sequence2, sequence.paths()[3:7])
        assert_is_iterable(sequence2)

        with pytest.raises(IndexError):
            _ = sequence[3:7:2]

    @staticmethod
    def tst_paths(sequence, filenames1):
        filenames2 = sequence.paths()
        for f1, f2 in zip(filenames1, filenames2):
            assert f1 == f2

    def tst_get_models(self, sequence, indices, outside_index):
        self.tst_get_models_index(sequence)
        for i in indices:
            self.tst_get_models_index(sequence, i)

    @staticmethod
    def tst_get_models_index(sequence, index=None):
        if index is None:
            sequence.get_detector()
            sequence.get_beam()
            sequence.get_goniometer()
            sequence.get_scan()
        else:
            sequence.get_detector(index)
            sequence.get_beam(index)
            sequence.get_goniometer(index)
            sequence.get_scan(index)

        # Ensure state at zero
        sequence[0]
        scan1 = sequence.get_scan()
        # Put sequence to end
        sequence[len(sequence) - 1]
        scan2 = sequence.get_scan()
        assert scan1 == scan2

    @staticmethod
    def tst_set_models(sequence):
        # Get some models
        beam = sequence.get_beam()
        gonio = sequence.get_goniometer()
        detector = sequence.get_detector()

        # Modify the geometry
        assert len(detector) == 1
        beam.set_direction((1, 0, 0))
        gonio.set_rotation_axis((0, 1, 0))
        detector[0].set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))

        # Override sequence models
        sequence.set_beam(beam)
        sequence.set_goniometer(gonio)
        sequence.set_detector(detector)

        # Ensure this doens't interfere with reading
        for i in sequence:
            pass

        # Get the models back and check they're ok
        beam2 = sequence.get_beam()
        gonio2 = sequence.get_goniometer()
        detector2 = sequence.get_detector()
        assert beam2 == beam
        assert gonio2 == gonio
        assert detector2 == detector

        # Get the models from an index back and check they're the same
        beam2 = sequence.get_beam(0)
        gonio2 = sequence.get_goniometer(0)
        detector2 = sequence.get_detector(0)
        assert beam2 == beam
        assert gonio2 == gonio
        assert detector2 == detector

        # Get a sub sequence
        sub_sequence = sequence[3:7]

        # Get the models back and check they're ok
        beam2 = sub_sequence.get_beam()
        gonio2 = sub_sequence.get_goniometer()
        detector2 = sub_sequence.get_detector()
        assert beam2 == beam
        assert gonio2 == gonio
        assert detector2 == detector

        # Get the models from an index back and check they're not the same
        beam2 = sub_sequence.get_beam(0)
        gonio2 = sub_sequence.get_goniometer(0)
        detector2 = sub_sequence.get_detector(0)
        assert beam2 == beam
        assert gonio2 == gonio
        assert detector2 == detector


@pytest.mark.parametrize("lazy", (True, False))
def test_SACLA_MPCCD_Cheetah_File(dials_regression, lazy):
    pytest.importorskip("h5py")
    filename = os.path.join(
        dials_regression,
        "image_examples",
        "SACLA_MPCCD_Cheetah",
        "run266702-0-subset.h5",
    )

    format_class = dxtbx.format.Registry.get_format_class_for_file(filename)

    iset = format_class.get_imageset([filename], lazy=lazy)

    assert len(iset) == 4
    for i in range(len(iset)):
        assert iset.get_raw_data(i)
        #   assert iset.get_mask(i)
        assert iset.get_beam(i)
        assert iset.get_detector(i)
        assert iset.get_goniometer(i) is None
        assert iset.get_scan(i) is None

    iset = format_class.get_imageset([filename], single_file_indices=[1], lazy=lazy)
    assert len(iset) == 1

    for i in range(len(iset)):
        assert iset.get_raw_data(i)
        #   assert iset.get_mask(i)
        assert iset.get_beam(i)
        assert iset.get_detector(i)
        assert iset.get_goniometer(i) is None
        assert iset.get_scan(i) is None


def test_imagesetfactory(centroid_files, dials_data):
    sequence = ImageSetFactory.new(centroid_files)

    assert isinstance(sequence[0], ImageSequence)

    template = dials_data("centroid_test_data").join("centroid_####.cbf").strpath
    image_range = (3, 6)

    sequence = ImageSetFactory.from_template(template, image_range)

    assert isinstance(sequence[0], ImageSequence)
    assert len(sequence[0]) == 4
    assert sequence[0].paths()[0].endswith("3.cbf")
    assert sequence[0].paths()[-1].endswith("6.cbf")

    imageset = ImageSetFactory.make_imageset(centroid_files)
    assert len(imageset) == 9

    imageset = ImageSetFactory.make_imageset(centroid_files, check_format=False)
    assert len(imageset) == 9

    sequence = ImageSetFactory.make_sequence(template, list(range(1, 9 + 1)))
    assert len(sequence) == 9

    sequence = ImageSetFactory.make_sequence(template, list(range(3, 6 + 1)))
    assert len(sequence) == 4


def test_pickle_imageset(centroid_files):
    sequence = ImageSetFactory.new(centroid_files)[0]

    # Read the 5th image
    assert sequence[4]

    # Pickle, then unpickle
    pickled_sequence = pickle.dumps(sequence)
    sequence2 = pickle.loads(pickled_sequence)

    assert sequence.get_template() == sequence2.get_template()
    assert sequence.get_array_range() == sequence2.get_array_range()
    assert sequence.get_beam() == sequence2.get_beam()
    assert sequence.get_goniometer() == sequence2.get_goniometer()
    assert sequence.get_scan() == sequence2.get_scan()
    assert sequence.paths() == sequence2.paths()
    assert sequence == sequence2

    # Check auxiliary methods after pickling
    sequence3 = sequence2[0:4]
    sequence4 = sequence3[0:2]
    sequence4.get_detectorbase(0)
    sequence4[0]
