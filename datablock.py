import collections
import errno
import itertools
import json
import logging
import math
import operator
import os.path
import pickle
from urllib.parse import urlparse

import libtbx
from scitbx import matrix

import dxtbx.imageset
import dxtbx.model
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.format.image import ImageBool, ImageDouble
from dxtbx.format.Registry import get_format_class_for_file
from dxtbx.sequence_filenames import (
    locate_files_matching_template_string,
    template_regex,
    template_string_number_index,
)
from dxtbx.serialize import load
from dxtbx.serialize.filename import resolve_path
from dxtbx.serialize.load import _decode_dict

try:
    from typing import Any, Callable, Dict, Generator, Iterable, List, Type
except ImportError:
    pass

logger = logging.getLogger(__name__)


def _iterate_with_previous(iterable):
    """Convenience iterator to give pairs of (previous, next) items"""
    previous = None
    for val in iterable:
        yield (previous, val)
        previous = val


class DataBlock:
    """High level container for blocks of sequences and imagesets."""

    def __init__(self, imagesets=None):
        """Instantiate from a list of imagesets."""
        # Try to get a format class
        self._format_class = None
        self._imagesets = []

        # Load imagesets
        if imagesets is not None:
            for iset in imagesets:
                self.append(iset)

    def append(self, imageset):
        """Add an imageset to the block."""
        if self._format_class is None:
            self._format_class = imageset.get_format_class()
        elif not self._format_class == imageset.get_format_class():
            raise TypeError("Can not mix image format classes in one datablock")
        self._imagesets.append(imageset)

    def extend(self, datablock):
        """Add two datablocks."""
        for iset in datablock:
            self.append(iset)

    def format_class(self):
        """Return the format class."""
        return self._format_class

    def extract_stills(self):
        """Extract all the still imagesets"""
        return list(self.iter_stills())

    def extract_sequences(self):
        """Extract all the sequences from the block."""
        return list(self.iter_sequences())

    def extract_imagesets(self):
        """Extract all imagesets."""
        return list(self._imagesets)

    def num_images(self):
        """Get the number of images."""
        return sum(len(iset) for iset in self._imagesets)

    def __len__(self):
        """The number of image sets."""
        return len(self._imagesets)

    def __eq__(self, rhs):
        """Check if two blocks are the same."""
        return (
            self._format_class == rhs._format_class
            and self._imagesets == rhs._imagesets
        )

    __hash__ = None  # datablock objects are mutable and therefore unhashable

    def __ne__(self, rhs):
        """Check if two blocks are not equal."""
        return not self.__eq__(rhs)

    def __iter__(self):
        """Iterate through the imagesets."""
        yield from self._imagesets

    def iter_sequences(self):
        """Iterate over sequence groups."""
        for iset in self._imagesets:
            if isinstance(iset, dxtbx.imageset.ImageSequence):
                yield iset

    def iter_stills(self):
        """Iterate over still groups."""
        for iset in self._imagesets:
            if not isinstance(iset, dxtbx.imageset.ImageSequence):
                yield iset

    def _find_unique_items(self, item_name, filter_none=False):
        """Return a list of unique beams, detectors, ... in order.
        Optionally filter out None values (unless they came in via
        an ImageSequence)."""

        # Use the keys of an ordered dictionary to guaranteeing uniqueness
        # while keeping the order. Could rewrite to use OrderedSet or wait
        # for Python 3.7 and use a regular dictionary
        # (3.7+: all dictionaries are guaranteed to be ordered)
        items = collections.OrderedDict()
        for imageset in self._imagesets:
            getter_function = getattr(imageset, "get_" + item_name)
            if isinstance(imageset, dxtbx.imageset.ImageSequence):
                items[getter_function()] = None
            else:
                for i in range(len(imageset)):
                    item = getter_function(i)
                    if not filter_none or item is not None:
                        items[item] = None
        return list(items)

    def unique_beams(self):
        return self._find_unique_items("beam")

    def unique_detectors(self):
        return self._find_unique_items("detector")

    def unique_goniometers(self):
        return self._find_unique_items("goniometer", filter_none=True)

    def unique_scans(self):
        return self._find_unique_items("scan", filter_none=True)

    def to_dict(self):
        """Convert the datablock to a dictionary"""

        def abspath_or_none(filename):
            if filename is None or filename == "":
                return None
            return os.path.abspath(filename)

        # Get a list of all the unique models
        b = list(self.unique_beams())
        d = list(self.unique_detectors())
        g = list(self.unique_goniometers())
        s = list(self.unique_scans())

        # Create the data block dictionary
        result = collections.OrderedDict()
        result["__id__"] = "DataBlock"
        result["imageset"] = []

        # Loop through all the imagesets
        for iset in self._imagesets:
            if isinstance(iset, dxtbx.imageset.ImageSequence):
                if iset.reader().is_single_file_reader():
                    result["imageset"].append(
                        collections.OrderedDict(
                            [
                                ("__id__", "ImageSequence"),
                                (
                                    "master",
                                    os.path.abspath(iset.reader().master_path()),
                                ),
                                (
                                    "mask",
                                    abspath_or_none(iset.external_lookup.mask.filename),
                                ),
                                (
                                    "gain",
                                    abspath_or_none(iset.external_lookup.gain.filename),
                                ),
                                (
                                    "pedestal",
                                    abspath_or_none(
                                        iset.external_lookup.pedestal.filename
                                    ),
                                ),
                                (
                                    "dx",
                                    abspath_or_none(iset.external_lookup.dx.filename),
                                ),
                                (
                                    "dy",
                                    abspath_or_none(iset.external_lookup.dy.filename),
                                ),
                                ("beam", b.index(iset.get_beam())),
                                ("detector", d.index(iset.get_detector())),
                                ("goniometer", g.index(iset.get_goniometer())),
                                ("scan", s.index(iset.get_scan())),
                                ("images", list(iset.indices())),
                                ("params", iset.params()),
                            ]
                        )
                    )
                else:
                    result["imageset"].append(
                        collections.OrderedDict(
                            [
                                ("__id__", "ImageSequence"),
                                ("template", os.path.abspath(iset.get_template())),
                                (
                                    "mask",
                                    abspath_or_none(iset.external_lookup.mask.filename),
                                ),
                                (
                                    "gain",
                                    abspath_or_none(iset.external_lookup.gain.filename),
                                ),
                                (
                                    "pedestal",
                                    abspath_or_none(
                                        iset.external_lookup.pedestal.filename
                                    ),
                                ),
                                (
                                    "dx",
                                    abspath_or_none(iset.external_lookup.dx.filename),
                                ),
                                (
                                    "dy",
                                    abspath_or_none(iset.external_lookup.dy.filename),
                                ),
                                ("beam", b.index(iset.get_beam())),
                                ("detector", d.index(iset.get_detector())),
                                ("goniometer", g.index(iset.get_goniometer())),
                                ("scan", s.index(iset.get_scan())),
                                ("params", iset.params()),
                            ]
                        )
                    )
            else:
                imageset = collections.OrderedDict()
                if isinstance(iset, dxtbx.imageset.ImageGrid):
                    imageset["__id__"] = "ImageGrid"
                    imageset["grid_size"] = iset.get_grid_size()
                else:
                    imageset["__id__"] = "ImageSet"
                image_list = []
                for i in range(len(iset)):
                    image_dict = collections.OrderedDict()
                    image_dict["filename"] = os.path.abspath(iset.get_path(i))
                    image_dict["gain"] = abspath_or_none(
                        iset.external_lookup.gain.filename
                    )
                    image_dict["pedestal"] = abspath_or_none(
                        iset.external_lookup.pedestal.filename
                    )
                    image_dict["dx"] = abspath_or_none(iset.external_lookup.dx.filename)
                    image_dict["dy"] = abspath_or_none(iset.external_lookup.dy.filename)
                    if iset.reader().is_single_file_reader():
                        image_dict["image"] = iset.indices()[i]
                    try:
                        image_dict["beam"] = b.index(iset.get_beam(i))
                    except Exception:
                        pass
                    try:
                        image_dict["detector"] = d.index(iset.get_detector())
                    except Exception:
                        pass
                    try:
                        image_dict["goniometer"] = g.index(iset.get_goniometer())
                    except Exception:
                        pass
                    try:
                        image_dict["scan"] = s.index(iset.get_scan(i))
                    except Exception:
                        pass
                    image_list.append(image_dict)
                imageset["mask"] = abspath_or_none(iset.external_lookup.mask.filename)
                imageset["images"] = image_list
                imageset["params"] = iset.params()
                result["imageset"].append(imageset)

        # Add the models to the dictionary
        result["beam"] = [bb.to_dict() for bb in b]
        result["detector"] = [dd.to_dict() for dd in d]
        result["goniometer"] = [gg.to_dict() for gg in g]
        result["scan"] = [ss.to_dict() for ss in s]

        return result


class FormatChecker:
    """A helper class to speed up identifying the correct image format by first
    trying the last format that was used."""

    def __init__(self):
        """Set the format class to none."""
        self._format_class = None

    def find_format(self, filename):
        """Search the registry for the image format class.
        Where possible use the last seen format class as a prioritisation hint.
        """
        if self._format_class:
            self._format_class = get_format_class_for_file(
                filename, format_hint=self._format_class.__name__
            )
        else:
            self._format_class = get_format_class_for_file(filename)
        if self._format_class:
            logger.debug("Using %s for %s", self._format_class.__name__, filename)
        else:
            logger.debug("No format class found for %s", filename)
        return self._format_class

    def iter_groups(self, filenames):
        group_format = None
        group_fnames = []
        for filename in filenames:
            fmt = self.find_format(filename)
            if fmt == group_format:
                group_fnames.append(filename)
            else:
                if group_fnames:
                    yield group_format, group_fnames
                group_fnames = [filename]
                group_format = fmt
            if fmt is not None:
                logger.debug("Using %s for %s", fmt.__name__, filename)
        if group_fnames:
            yield group_format, group_fnames


class DataBlockTemplateImporter:
    """A class to import a datablock from a template."""

    def __init__(self, templates, **kwargs):
        """Import the datablocks from the given templates."""
        assert "verbose" not in kwargs, "The verbose parameter has been removed"
        assert len(templates) > 0

        self.datablocks = []

        # A function to append or create a new datablock
        def append_to_datablocks(iset):
            try:
                self.datablocks[-1].append(iset)
            except (IndexError, TypeError):
                # This happens when we already have a datablock with a different format
                self.datablocks.append(DataBlock([iset]))
            logger.debug("Added imageset to datablock %d", len(self.datablocks) - 1)

        # For each template do an import
        for template in templates:
            template = os.path.normpath(template)
            paths = sorted(locate_files_matching_template_string(template))
            logger.debug("The following files matched the template string:")
            if len(paths) > 0:
                for p in paths:
                    logger.debug(" %s", p)
            else:
                logger.debug(" No files found")

            # Check if we've matched any filenames
            if len(paths) == 0:
                raise ValueError('Template "%s" does not match any files' % template)

            # Get the format from the first image
            fmt = FormatChecker().find_format(paths[0])
            if fmt is None:
                raise ValueError("Image file %s format is unknown" % paths[0])
            elif fmt.is_abstract():
                raise ValueError(
                    f"Image file {paths[0]} appears to be a '{type(fmt).__name__}', but this is an abstract Format"
                )
            else:
                imageset = self._create_imageset(fmt, template, paths, **kwargs)
                append_to_datablocks(imageset)

    def _create_imageset(self, format_class, template, filenames, **kwargs):
        """Create a multi file sequence or imageset."""
        # Get the image range
        index = slice(*template_string_number_index(template))
        image_range = kwargs.get("image_range")
        if image_range:
            first, last = image_range
        else:
            first = int(filenames[0][index])
            last = int(filenames[-1][index])

        # Check all images in range are present - if allowed
        if not kwargs.get("allow_incomplete_sequences", False):
            all_numbers = {int(f[index]) for f in filenames}
            missing = set(range(first, last + 1)) - all_numbers
            if missing:
                raise ValueError(
                    "Missing image{} {} from imageset ({}-{})".format(
                        "s" if len(missing) > 1 else "",
                        ", ".join(str(x) for x in sorted(missing)),
                        first,
                        last,
                    )
                )

        # Read the image
        fmt = format_class(filenames[0], **(kwargs.get("format_kwargs", {})))

        # Get the meta data from the format
        b = fmt.get_beam()
        d = fmt.get_detector()
        g = fmt.get_goniometer()
        s = fmt.get_scan()

        # Update the image range
        s.set_image_range((first, last))

        # Get the image range
        image_range = s.get_image_range()
        image_range = (image_range[0], image_range[1] + 1)

        # Create the sequence
        imageset = dxtbx.imageset.ImageSetFactory.make_sequence(
            template,
            list(range(*image_range)),
            format_class,
            b,
            d,
            g,
            s,
            format_kwargs=kwargs.get("format_kwargs"),
        )

        return imageset


class ImageMetadataRecord:
    """Object to store metadata information.

    This is used whilst building the datablocks.  The metadata for each
    image can be read once, and then any grouping/deduplication can happen
    later, without re-opening the original file.
    """

    def __init__(
        self,
        beam=None,
        detector=None,
        goniometer=None,
        scan=None,
        template=None,
        filename=None,
        index=None,
    ):
        # type: (dxtbx.model.Beam, dxtbx.model.Detector, dxtbx.model.Goniometer, dxtbx.model.Scan, str, str, int)
        """
        Args:
            beam:       Stores a beam model
            detector:   Stores a detector model
            goniometer: Stores a goniometer model
            scan:       Stores a scan model
            filename:   The filename this record was parsed from
            template:
                The template string parsed from the filename. Usually,
                the template is only present if a scan was found and
                oscillation width was nonzero.
            index:
                The index of this file in the template. Applying the
                index to the template field should recover the filename
        """
        self.beam = beam
        self.detector = detector
        self.goniometer = goniometer
        self.scan = scan
        self.template = template
        self.filename = filename
        self.index = index

    def merge_metadata_from(
        self,
        other_record,
        compare_beam=operator.__eq__,
        compare_detector=operator.__eq__,
        compare_goniometer=operator.__eq__,
    ):
        # type: (ImageMetadataRecord, Callable, Callable, Callable) -> bool
        """
        Compare two record objects and merge equivalent data.

        This method will compare (with optional functions) instance data
        for beam, detector and goniometer. If any of the metadata for
        this record is equivalent to (but a different instance from) the
        other record, then this instance will be altered to match the
        other. The function used to compare beams, detectors and
        goniometers can be customised - but by default the normal
        equality operator is used.

        Args:
            other_record:       Another metadata instance
            compare_beam:       A function to compare beams
            compare_detector:   A function to compare detectors
            compare_goniometer: A function to compare goniometers

        Returns: True if any action was taken
        """
        # Allow 'defaults' of None to work - behavior from legacy implementation
        compare_beam = compare_beam or operator.__eq__
        compare_detector = compare_detector or operator.__eq__
        compare_goniometer = compare_goniometer or operator.__eq__

        record_altered = False
        if self.beam is not other_record.beam and compare_beam(
            self.beam, other_record.beam
        ):
            self.beam = other_record.beam
            record_altered = True
        if self.detector is not other_record.detector and compare_detector(
            self.detector, other_record.detector
        ):
            self.detector = other_record.detector
            record_altered = True
        if self.goniometer is not other_record.goniometer and compare_goniometer(
            self.goniometer, other_record.goniometer
        ):
            self.goniometer = other_record.goniometer
            record_altered = True

        return record_altered

    @classmethod
    def from_format(cls, fmt):
        # type: (Format) -> Any
        """
        Read metadata information from a Format instance.

        This will only pull information out of a single format instance
        while it is open - combining metadata records must be done
        separately.

        Args:
            format: The instance of the format class to read data from

        Returns:
            A new ImageMetadataRecord with the pre-read information
        """
        record = cls()
        record.filename = fmt.get_image_file()
        # Get the metadata from the format
        try:
            record.beam = fmt.get_beam()
        except Exception:
            pass
        try:
            record.detector = fmt.get_detector()
        except Exception:
            pass
        try:
            record.goniometer = fmt.get_goniometer()
        except Exception:
            pass
        try:
            record.scan = fmt.get_scan()
        except Exception:
            pass

        # Get the template and index if possible - and only if we've got a
        # recorded oscillation value
        if record.scan is not None:
            record.template, record.index = template_regex(record.filename)

        return record

    def __repr__(self):
        items = [
            ("filename", self.filename),
            ("beam", self.beam),
            ("detector", self.detector),
            ("goiometer", self.goniometer),
            ("scan", self.scan),
            ("template", self.template),
            ("index", self.index),
        ]
        itemstr = ", ".join(x + "=" + repr(y) for x, y in items)
        return "<{}{}{}>".format(type(self).__name__, " " if itemstr else "", itemstr)

    def __hash__(self):
        return hash(
            (
                self.beam,
                self.detector,
                self.goniometer,
                self.scan,
                self.template,
                self.filename,
                self.index,
            )
        )

    def __eq__(self, other):
        if not isinstance(other, ImageMetadataRecord):
            return False
        return all(
            getattr(self, attribute) == getattr(other, attribute)
            for attribute in (
                "beam",
                "detector",
                "goniometer",
                "scan",
                "template",
                "filename",
                "index",
            )
        )

    def __ne__(self, other):
        return not self == other


def _openingpathiterator(pathnames: Iterable[str]):
    """Utility function to efficiently open all paths.

    A path is a potential file or directory.
    Each path will be opened with :meth:`dxtbx.format.Format.open_file`,
    but in order to do so each file will only be opened once, and extraneous
    use of :func:`os.stat` will be avoided.
    Any path entries that are a directory will be recursed into, once -
    any further directories found will be ignored. Any path that is not
    a file or directory, or on which IO fails for any reason, will still
    be returned.

    Args:
        pathnames: Paths to attempt to open
    """

    # Store a tuple of (recurse, pathname) to track what was root level
    paths = collections.deque((True, x) for x in sorted(pathnames))

    while paths:
        # Get the next path from the queue
        (do_recurse, pathname) = paths.popleft()

        try:
            # Attempt to open this 'path'
            Format.open_file(pathname)
        except OSError as e:
            if e.errno == errno.EISDIR:
                if do_recurse:
                    # We've tried to open a directory. Get all the entries...
                    subdir_paths = sorted(
                        os.path.join(pathname, x) for x in os.listdir(pathname)
                    )
                    # ... and add them to our queue. Make sure not to mark for recursion
                    paths.extendleft((False, x) for x in reversed(subdir_paths))
                    logger.debug("Adding %d files from %s", len(subdir_paths), pathname)
                else:
                    logger.debug("Not adding sub-level directory entry %s", pathname)
                # Don't return directory instances
                continue
            else:
                # A non-directory-related IO error
                logger.debug("Could not import %s: %s", pathname, os.strerror(e.errno))

        yield pathname


def _merge_model_metadata(
    records, compare_beam=None, compare_detector=None, compare_goniometer=None
):
    # type: (Iterable[ImageMetadataRecord], Callable, Callable, Callable)
    """
    Merge metadata between consecutive record objects.

    This will compare each record with the previous one, and make sure
    the metadata instances are shared where appropriate.

    Args:
        records:    Records for the images to merge into imagesets
        compare_beam:       The function to to compare beams
        compare_detector:   The function to compare detectors
        compare_goniometer: The function to compare goniometers
    """
    for prev, record in _iterate_with_previous(records):
        if prev is None:
            continue
        record.merge_metadata_from(
            prev,
            compare_beam=compare_beam,
            compare_detector=compare_detector,
            compare_goniometer=compare_goniometer,
        )


def _merge_scans(records, scan_tolerance=None):
    # type: (Iterable[ImageMetadataRecord], float) -> List[ImageMetadataRecord]
    """
    Merge consecutive scan records with identical metadata.

    The records should have previously had their model metadata merged,
    as identity will be used to compare metadata identity at this stage.

    Args:
        records:        Records to merge
        scan_tolerance: Percentage of oscillation range to tolerate
                        when merging scan records

    Returns:
        A (potentially shorter) list of records with scans merged
    """
    merged_records = []
    logger.debug("Merging scans")
    for prev, record in _iterate_with_previous(records):
        # The first record always gets recorded
        if prev is None:
            merged_records.append(record)
            logger.debug("  Saving initial record %s", record)
            continue
        # Compare metadata instances
        same_metadata = [
            prev.beam is record.beam,
            prev.detector is record.detector,
            prev.goniometer is record.goniometer,
        ]
        # Condition for combining:
        # - All metadata must match
        # - Previous record must be templated
        # - This record must be templated
        if (
            all(same_metadata)
            and prev.template is not None
            and record.template is not None
        ):
            # Attempt to append to scan
            try:
                if scan_tolerance is None:
                    prev.scan.append(record.scan)
                else:
                    prev.scan.append(record.scan, scan_tolerance=scan_tolerance)
            except RuntimeError as e:
                print(e)
                logger.debug(
                    "  Failed to merge record %s with previous - writing new scan"
                )
            else:
                # If we appended, then we don't need to keep this record's scan
                record.scan = prev.scan
                logger.debug("  Appended record %s to previous", record)
                continue
        merged_records.append(record)

    logger.debug("Result of merging record scans: %d records", len(merged_records))
    return merged_records


def _create_imagesequence(record, format_class, format_kwargs=None):
    # type: (ImageMetadataRecord, Type[Format], Dict) -> dxtbx.imageset.ImageSequence
    """
    Create an ImageSequence object from a single rotation data image.

    Args:
        record: Single-image metadata records to merge into a single imageset
        format_class: The format class object for these image records
        format_kwargs: Extra arguments to pass to the format class when
            creating an ImageSet

    Returns:
        An imageset representing the sequence of data
    """
    index_start, index_end = record.scan.get_image_range()
    # Create the sequence
    sequence = dxtbx.imageset.ImageSetFactory.make_sequence(
        template=os.path.abspath(record.template),
        indices=list(range(index_start, index_end + 1)),
        format_class=format_class,
        beam=record.beam,
        detector=record.detector,
        goniometer=record.goniometer,
        scan=record.scan,
        format_kwargs=format_kwargs,
        # check_format=False,
    )
    return sequence


def _groupby_template_is_none(records):
    # type: (Iterable[ImageMetadataRecord]) -> Generator[List[ImageMetadataRecord]]
    """Specialization of groupby that groups records by format=None"""
    for _, group in itertools.groupby(
        enumerate(records), key=lambda x: -1 if x[1].template is None else x[0]
    ):
        yield list(x[1] for x in group)


def _convert_to_imagesets(records, format_class, format_kwargs=None):
    # type: (Iterable[ImageMetadataRecord], Type[dxtbx.format.Format], Dict) -> Generator[dxtbx.imageset.ImageSet]
    """
    Convert records into imagesets.

    The records should have been metadata- and scan-merged by this point.
    Rules:
    - Any groups of template=None where any of the metadata objects
      are shared, go into a single imageset
    - Anything with a template goes into a single sequence

    Args:
        records: The records to convert
        format_class: The format class for the data in this record
        format_kwargs: Any format configuration arguments to pass
            to the format imageset creator.

    Returns:
        Imagesets representing the records
    """

    # Iterate over images/sets such that template=None are clustered
    for setgroup in _groupby_template_is_none(records):
        if setgroup[0].template is not None:
            # If we have a template, then it's a sequence
            assert len(setgroup) == 1, "Got group of metadata records in template?"
            logger.debug("Creating Imagesequence from %s", setgroup[0].template)
            yield _create_imagesequence(setgroup[0], format_class, format_kwargs)
        else:
            # Without a template, it was never identified as a sequence, so an imageset
            logger.debug("Creating ImageSet from %d files", len(setgroup))
            yield _create_imageset(setgroup, format_class, format_kwargs)


def _create_imageset(records, format_class, format_kwargs=None):
    # type: (Iterable[ImageMetadataRecord], Type[dxtbx.format.Format], Dict) -> dxtbx.imageset.ImageSet
    """
    Create an ImageSet object from a set of single-image records.

    Args:
        records: Single-image metadata records to merge into a single imageset
        format_class: The format class object for these image records
        format_kwargs: Extra arguments to pass to the format class when
            creating an ImageSet

    Returns:
        An imageset for all the image records
    """
    records = list(records)
    # Nothing here should have been assigned a template parameter
    assert all(x.template is None for x in records)
    # Extract the filenames from the records
    filenames = [
        x.filename if urlparse(x.filename).scheme else os.path.abspath(x.filename)
        for x in records
    ]
    # Create the imageset
    imageset = dxtbx.imageset.ImageSetFactory.make_imageset(
        filenames, format_class, format_kwargs=format_kwargs, check_format=False
    )
    # Update all of the metadata for each record
    for i, r in enumerate(records):
        imageset.set_beam(r.beam, i)
        imageset.set_detector(r.detector, i)
        imageset.set_goniometer(r.goniometer, i)
        imageset.set_scan(r.scan, i)
    return imageset


class DataBlockFilenameImporter:
    """A class to import a datablock from image files."""

    def __init__(
        self,
        filenames,
        compare_beam=None,
        compare_detector=None,
        compare_goniometer=None,
        scan_tolerance=None,
        format_kwargs=None,
    ):
        """Import the datablocks from the given filenames."""

        # Init the datablock list
        self.unhandled = []
        self.datablocks = []

        # A function to append or create a new datablock
        def append_to_datablocks(iset):
            try:
                self.datablocks[-1].append(iset)
            except (IndexError, TypeError):
                # This happens when we already have a datablock with a different format
                self.datablocks.append(DataBlock([iset]))
            logger.debug("Added imageset to datablock %d", len(self.datablocks) - 1)

        # Process each file given by this path list
        to_process = _openingpathiterator(filenames)
        find_format = FormatChecker()

        format_groups = collections.OrderedDict()
        if format_kwargs is None:
            format_kwargs = {}
        for filename in to_process:
            # We now have a file, pre-opened by Format.open_file (therefore
            # cached). Determine its type, and prepare to put into a group
            format_class = find_format.find_format(filename)

            # Verify this makes sense
            if not format_class:
                # No format class found?
                logger.debug("Could not determine format for %s", filename)
                self.unhandled.append(filename)
            elif format_class.is_abstract():
                logger.debug(
                    f"Image file {filename} appears to be a '{format_class.__name__}', but this is an abstract Format"
                )
                # Invalid format class found?
                self.unhandled.append(filename)
            elif issubclass(format_class, FormatMultiImage):
                imageset = self._create_single_file_imageset(
                    format_class, filename, format_kwargs=format_kwargs
                )
                format_groups.setdefault(format_class, []).append(imageset)
                logger.debug("Loaded file: %s", filename)
            else:
                format_object = format_class(filename, **format_kwargs)
                meta = ImageMetadataRecord.from_format(format_object)
                assert meta.filename == filename

                # Add this entry to our table of formats
                format_groups.setdefault(format_class, []).append(meta)
                logger.debug("Loaded metadata of file: %s", filename)

        # Now, build datablocks from these files. Duplicating the logic of
        # the previous implementation:
        # - Each change in format goes into its own Datablock
        # - FormatMultiImage files each have their own ImageSet
        # - Every set of images forming a scan goes into its own ImageSequence
        # - Any consecutive still frames that share any metadata with the
        #   previous still fram get collected into one ImageSet

        # Treat each format as a separate datablock
        for format_class, records in format_groups.items():
            if issubclass(format_class, FormatMultiImage):
                for imageset in records:
                    append_to_datablocks(imageset)
                continue

            # Merge any consecutive and identical metadata together
            _merge_model_metadata(
                records,
                compare_beam=compare_beam,
                compare_detector=compare_detector,
                compare_goniometer=compare_goniometer,
            )
            records = _merge_scans(records, scan_tolerance=scan_tolerance)
            imagesets = _convert_to_imagesets(records, format_class, format_kwargs)
            imagesets = list(imagesets)

            # Validate this datablock and store it
            assert imagesets, "Datablock got no imagesets?"
            for i in imagesets:
                append_to_datablocks(i)

        # Now we are done. since in an __init__, the caller can now pick the
        # results out of the .datablocks and .unhandled instance attributes.

    def _create_single_file_imageset(self, format_class, filename, format_kwargs=None):
        """Create an imageset from a multi image file."""
        if format_kwargs is None:
            format_kwargs = {}
        return format_class.get_imageset(
            os.path.abspath(filename), format_kwargs=format_kwargs
        )


class InvalidDataBlockError(RuntimeError):
    """
    Indicates an error whilst validating the experiment list.

    This means that there is some structural problem that prevents the given data
    from representing a well-formed experiment list. This doesn't indicate e.g.
    some problem with the data or model consistency.
    """


class DataBlockDictImporter:
    """A class to import a datablock from dictionary."""

    def __init__(self, obj, check_format=True, directory=None):
        """Get the datablocks from the dictionary."""
        self.datablocks = self._load_datablocks(obj, check_format, directory)

    def _load_datablocks(self, obj, check_format=True, directory=None):
        """Create the datablock from a dictionary."""

        # If we have a list, extract for each dictionary in the list
        if isinstance(obj, list):
            return [self._load_datablocks(dd, check_format, directory) for dd in obj]
        elif not isinstance(obj, dict):
            raise InvalidDataBlockError(
                "Unexpected datablock type {} instead of dict".format(type(obj))
            )
        # Make sure the id signature is correct
        if not obj.get("__id__") == "DataBlock":
            raise InvalidDataBlockError(
                "Expected __id__ 'DataBlock', but found {}".format(
                    repr(obj.get("__id__"))
                )
            )

        # Get the list of models
        blist = obj.get("beam", [])
        dlist = obj.get("detector", [])
        glist = obj.get("goniometer", [])
        slist = obj.get("scan", [])

        def load_models(obj):
            try:
                beam = dxtbx.model.BeamFactory.from_dict(blist[obj["beam"]])
            except Exception:
                beam = None
            try:
                dobj = dlist[obj["detector"]]
                detector = dxtbx.model.DetectorFactory.from_dict(dobj)
            except Exception:
                detector = None
            try:
                gonio = dxtbx.model.GoniometerFactory.from_dict(
                    glist[obj["goniometer"]]
                )
            except Exception:
                gonio = None
            try:
                scan = dxtbx.model.ScanFactory.from_dict(slist[obj["scan"]])
            except Exception:
                scan = None
            return beam, detector, gonio, scan

        # Loop through all the imagesets
        imagesets = []
        for imageset in obj["imageset"]:
            ident = imageset["__id__"]
            if "params" in imageset:
                format_kwargs = imageset["params"]
            else:
                format_kwargs = {}
            if ident == "ImageSequence" or ident == "ImageSweep":
                beam, detector, gonio, scan = load_models(imageset)
                if "template" in imageset:
                    template = resolve_path(imageset["template"], directory=directory)
                    i0, i1 = scan.get_image_range()
                    iset = dxtbx.imageset.ImageSetFactory.make_sequence(
                        template,
                        list(range(i0, i1 + 1)),
                        None,
                        beam,
                        detector,
                        gonio,
                        scan,
                        check_format,
                        format_kwargs=format_kwargs,
                    )
                    if "mask" in imageset and imageset["mask"] is not None:
                        imageset["mask"] = resolve_path(
                            imageset["mask"], directory=directory
                        )
                        iset.external_lookup.mask.filename = imageset["mask"]
                        if check_format:
                            with open(imageset["mask"], "rb") as infile:
                                iset.external_lookup.mask.data = ImageBool(
                                    pickle.load(infile, encoding="bytes")
                                )
                    if "gain" in imageset and imageset["gain"] is not None:
                        imageset["gain"] = resolve_path(
                            imageset["gain"], directory=directory
                        )
                        iset.external_lookup.gain.filename = imageset["gain"]
                        if check_format:
                            with open(imageset["gain"], "rb") as infile:
                                iset.external_lookup.gain.data = ImageDouble(
                                    pickle.load(infile, encoding="bytes")
                                )
                    if "pedestal" in imageset and imageset["pedestal"] is not None:
                        imageset["pedestal"] = resolve_path(
                            imageset["pedestal"], directory=directory
                        )
                        iset.external_lookup.pedestal.filename = imageset["pedestal"]
                        if check_format:
                            with open(imageset["pedestal"], "rb") as infile:
                                iset.external_lookup.pedestal.data = ImageDouble(
                                    pickle.load(infile, encoding="bytes")
                                )
                    if "dx" in imageset and imageset["dx"] is not None:
                        imageset["dx"] = resolve_path(
                            imageset["dx"], directory=directory
                        )
                        iset.external_lookup.dx.filename = imageset["dx"]
                        with open(imageset["dx"], "rb") as infile:
                            iset.external_lookup.dx.data = ImageDouble(
                                pickle.load(infile, encoding="bytes")
                            )
                    if "dy" in imageset and imageset["dy"] is not None:
                        imageset["dy"] = resolve_path(
                            imageset["dy"], directory=directory
                        )
                        iset.external_lookup.dy.filename = imageset["dy"]
                        with open(imageset["dy"], "rb") as infile:
                            iset.external_lookup.dy.data = ImageDouble(
                                pickle.load(infile, encoding="bytes")
                            )
                    iset.update_detector_px_mm_data()
                elif "master" in imageset:
                    template = resolve_path(imageset["master"], directory=directory)
                    i0, i1 = scan.get_image_range()
                    if not check_format:
                        format_class = FormatMultiImage
                    else:
                        format_class = None
                    iset = dxtbx.imageset.ImageSetFactory.make_sequence(
                        template,
                        list(range(i0, i1 + 1)),
                        format_class=format_class,
                        beam=beam,
                        detector=detector,
                        goniometer=gonio,
                        scan=scan,
                        check_format=check_format,
                        format_kwargs=format_kwargs,
                    )
                    if "mask" in imageset and imageset["mask"] is not None:
                        imageset["mask"] = resolve_path(imageset["mask"], directory)
                        iset.external_lookup.mask.filename = imageset["mask"]
                        if check_format:
                            with open(imageset["mask"], "rb") as infile:
                                iset.external_lookup.mask.data = ImageBool(
                                    pickle.load(infile, encoding="bytes")
                                )
                    if "gain" in imageset and imageset["gain"] is not None:
                        imageset["gain"] = resolve_path(imageset["gain"], directory)
                        iset.external_lookup.gain.filename = imageset["gain"]
                        if check_format:
                            with open(imageset["gain"], "rb") as infile:
                                iset.external_lookup.gain.data = ImageDouble(
                                    pickle.load(infile, encoding="bytes")
                                )
                    if "pedestal" in imageset and imageset["pedestal"] is not None:
                        imageset["pedestal"] = resolve_path(
                            imageset["pedestal"], directory
                        )
                        iset.external_lookup.pedestal.filename = imageset["pedestal"]
                        if check_format:
                            with open(imageset["pedestal"], "rb") as infile:
                                iset.external_lookup.pedestal.data = ImageDouble(
                                    pickle.load(infile, encoding="bytes")
                                )
                    if "dx" in imageset and imageset["dx"] is not None:
                        imageset["dx"] = resolve_path(imageset["dx"], directory)
                        iset.external_lookup.dx.filename = imageset["dx"]
                        with open(imageset["dx"], "rb") as infile:
                            iset.external_lookup.dx.data = ImageDouble(
                                pickle.load(infile, encoding="bytes")
                            )
                    if "dy" in imageset and imageset["dy"] is not None:
                        imageset["dy"] = resolve_path(imageset["dy"], directory)
                        iset.external_lookup.dy.filename = imageset["dy"]
                        with open(imageset["dy"], "rb") as infile:
                            iset.external_lookup.dy.data = ImageDouble(
                                pickle.load(infile, encoding="bytes")
                            )
                    iset.update_detector_px_mm_data()
                imagesets.append(iset)
            elif ident == "ImageSet" or ident == "ImageGrid":
                filenames = [image["filename"] for image in imageset["images"]]
                indices = [
                    image["image"] for image in imageset["images"] if "image" in image
                ]
                assert len(indices) == 0 or len(indices) == len(filenames)
                iset = dxtbx.imageset.ImageSetFactory.make_imageset(
                    filenames, None, check_format, indices, format_kwargs=format_kwargs
                )
                if ident == "ImageGrid":
                    grid_size = imageset["grid_size"]
                    iset = dxtbx.imageset.ImageGrid.from_imageset(iset, grid_size)
                for i, image in enumerate(imageset["images"]):
                    beam, detector, gonio, scan = load_models(image)
                    iset.set_beam(beam, i)
                    iset.set_detector(detector, i)
                    iset.set_goniometer(gonio, i)
                    iset.set_scan(scan, i)
                if "mask" in imageset and imageset["mask"] is not None:
                    imageset["mask"] = resolve_path(imageset["mask"], directory)
                    iset.external_lookup.mask.filename = imageset["mask"]
                    if check_format:
                        with open(imageset["mask"], "rb") as infile:
                            iset.external_lookup.mask.data = ImageBool(
                                pickle.load(infile, encoding="bytes")
                            )
                if "gain" in imageset and imageset["gain"] is not None:
                    imageset["gain"] = resolve_path(imageset["gain"], directory)
                    iset.external_lookup.gain.filename = imageset["gain"]
                    if check_format:
                        with open(imageset["gain"], "rb") as infile:
                            iset.external_lookup.gain.data = ImageDouble(
                                pickle.load(infile, encoding="bytes")
                            )
                if "pedestal" in imageset and imageset["pedestal"] is not None:
                    imageset["pedestal"] = resolve_path(imageset["pedestal"], directory)
                    iset.external_lookup.pedestal.filename = imageset["pedestal"]
                    if check_format:
                        with open(imageset["pedestal"], "rb") as infile:
                            iset.external_lookup.pedestal.data = ImageDouble(
                                pickle.load(infile, encoding="bytes")
                            )
                if "dx" in imageset and imageset["dx"] is not None:
                    imageset["dx"] = resolve_path(imageset["dx"], directory)
                    iset.external_lookup.dx.filename = imageset["dx"]
                    with open(imageset["dx"], "rb") as infile:
                        iset.external_lookup.dx.data = ImageDouble(
                            pickle.load(infile, encoding="bytes")
                        )
                if "dy" in imageset and imageset["dy"] is not None:
                    imageset["dy"] = resolve_path(imageset["dy"], directory)
                    iset.external_lookup.dy.filename = imageset["dy"]
                    with open(imageset["dy"], "rb") as infile:
                        iset.external_lookup.dy.data = ImageDouble(
                            pickle.load(infile, encoding="bytes")
                        )
                    iset.update_detector_px_mm_data()
                imagesets.append(iset)
            else:
                raise RuntimeError("expected ImageSet/ImageSequence, got %s" % ident)

        return DataBlock(imagesets)


class DataBlockImageSetImporter:
    """A class to import a datablock from imagesets."""

    def __init__(self, imagesets):
        """Load a list of datablocks from imagesets."""
        self.datablocks = []
        if not isinstance(imagesets, list):
            imagesets = [imagesets]
        for imageset in imagesets:
            try:
                self.datablocks[-1].append(imageset)
            except (IndexError, AssertionError):
                self.datablocks.append(DataBlock([imageset]))


class DataBlockFactory:
    """Class for creating DataBlock instances"""

    @staticmethod
    def from_args(
        args,
        verbose=None,
        unhandled=None,
        compare_beam=None,
        compare_detector=None,
        compare_goniometer=None,
        scan_tolerance=None,
        format_kwargs=None,
    ):
        """Try to load datablocks from any recognized format."""

        if unhandled is None:
            unhandled = []
        unhandled1 = []

        # Try as image files
        datablocks = DataBlockFactory.from_filenames(
            args,
            unhandled=unhandled1,
            compare_beam=compare_beam,
            compare_detector=compare_detector,
            compare_goniometer=compare_goniometer,
            scan_tolerance=scan_tolerance,
            format_kwargs=format_kwargs,
        )

        # Try as serialized formats
        for filename in unhandled1:
            try:
                datablocks.extend(DataBlockFactory.from_serialized_format(filename))
                logger.debug("Loaded datablocks(s) from %s", filename)
            except Exception:
                unhandled.append(filename)

        return datablocks

    @staticmethod
    def from_filenames(
        filenames,
        verbose=None,
        unhandled=None,
        compare_beam=None,
        compare_detector=None,
        compare_goniometer=None,
        scan_tolerance=None,
        format_kwargs=None,
    ):
        """Create a list of data blocks from a list of directory or file names."""
        importer = DataBlockFilenameImporter(
            filenames,
            compare_beam=compare_beam,
            compare_detector=compare_detector,
            compare_goniometer=compare_goniometer,
            scan_tolerance=scan_tolerance,
            format_kwargs=format_kwargs,
        )
        if unhandled is not None:
            unhandled.extend(importer.unhandled)
        return importer.datablocks

    @staticmethod
    def from_dict(obj, check_format=True, directory=None):
        """Create a datablock from a dictionary."""
        importer = DataBlockDictImporter(obj, check_format, directory)
        return importer.datablocks

    @staticmethod
    def from_json(string, check_format=True, directory=None):
        """Decode a datablock from JSON string."""
        return DataBlockFactory.from_dict(
            json.loads(string, object_hook=_decode_dict),
            check_format=check_format,
            directory=directory,
        )

    @staticmethod
    def from_json_file(filename, check_format=True):
        """Decode a datablock from a JSON file."""
        filename = os.path.abspath(filename)
        directory = os.path.dirname(filename)
        with open(filename) as infile:
            return DataBlockFactory.from_json(
                infile.read(), check_format=check_format, directory=directory
            )

    @staticmethod
    def from_pickle_file(filename):
        """Decode a datablock from a pickle file."""
        with open(filename, "rb") as infile:
            obj = pickle.load(infile)
            if isinstance(obj, list):
                assert all(isinstance(db, DataBlock) for db in obj)
            else:
                assert isinstance(obj, DataBlock)
            return obj

    @staticmethod
    def from_imageset(imagesets):
        """Load a datablock from a list of imagesets."""
        importer = DataBlockImageSetImporter(imagesets)
        return importer.datablocks

    @staticmethod
    def from_imageset_json_file(filename):
        """Load a datablock from a sequence file."""
        # Load the imageset and create a datablock from the filenames
        imageset = load.imageset(filename)
        return DataBlockFactory.from_imageset(imageset)

    @staticmethod
    def from_serialized_format(filename, check_format=True):
        """Load a datablock from serialized formats."""

        # First try as JSON format
        try:
            return DataBlockFactory.from_json_file(filename, check_format)
        except Exception:
            pass

        # Now try as pickle format
        try:
            return DataBlockFactory.from_pickle_file(filename)
        except Exception:
            pass

        # Now try as imageset json files
        return DataBlockFactory.from_imageset_json_file(filename)

    @staticmethod
    def from_in_memory(images, indices=None):
        """Function to instantiate data block from in memory imageset."""
        return DataBlock(
            [
                dxtbx.imageset.ImageSet(
                    dxtbx.imageset.ImageSetData(dxtbx.imageset.MemReader(images)),
                    indices,
                )
            ]
        )


class AutoEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, libtbx.AutoType):
            return "Auto"
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


class BeamComparison:
    """A class to provide simple beam comparison"""

    def __init__(
        self,
        wavelength_tolerance=1e-6,
        direction_tolerance=1e-6,
        polarization_normal_tolerance=1e-6,
        polarization_fraction_tolerance=1e-6,
    ):
        self.wavelength_tolerance = wavelength_tolerance
        self.direction_tolerance = direction_tolerance
        self.polarization_normal_tolerance = polarization_normal_tolerance
        self.polarization_fraction_tolerance = polarization_fraction_tolerance

    def __call__(self, a, b):
        if a is None and b is None:
            return True
        return a.is_similar_to(
            b,
            wavelength_tolerance=self.wavelength_tolerance,
            direction_tolerance=self.direction_tolerance,
            polarization_normal_tolerance=self.polarization_normal_tolerance,
            polarization_fraction_tolerance=self.polarization_fraction_tolerance,
        )


class DetectorComparison:
    """A class to provide simple detector comparison"""

    def __init__(
        self, fast_axis_tolerance=1e-6, slow_axis_tolerance=1e-6, origin_tolerance=1e-6
    ):
        self.fast_axis_tolerance = fast_axis_tolerance
        self.slow_axis_tolerance = slow_axis_tolerance
        self.origin_tolerance = origin_tolerance

    def __call__(self, a, b):
        if a is None and b is None:
            return True
        return a.is_similar_to(
            b,
            fast_axis_tolerance=self.fast_axis_tolerance,
            slow_axis_tolerance=self.slow_axis_tolerance,
            origin_tolerance=self.origin_tolerance,
        )


class GoniometerComparison:
    """A class to provide simple goniometer comparison"""

    def __init__(
        self,
        rotation_axis_tolerance=1e-6,
        fixed_rotation_tolerance=1e-6,
        setting_rotation_tolerance=1e-6,
    ):
        self.rotation_axis_tolerance = rotation_axis_tolerance
        self.fixed_rotation_tolerance = fixed_rotation_tolerance
        self.setting_rotation_tolerance = setting_rotation_tolerance

    def __call__(self, a, b):
        if a is None and b is None:
            return True
        elif a is None or b is None:
            return False
        return a.is_similar_to(
            b,
            rotation_axis_tolerance=self.rotation_axis_tolerance,
            fixed_rotation_tolerance=self.fixed_rotation_tolerance,
            setting_rotation_tolerance=self.setting_rotation_tolerance,
        )


def _all_equal(a, b):
    return all(x[0] == x[1] for x in zip(a, b))


def _all_approx_equal(a, b, tolerance):
    return all(abs(x[0] - x[1]) < tolerance for x in zip(a, b))


class BeamDiff:
    """A class to provide simple beam comparison"""

    def __init__(
        self,
        wavelength_tolerance=1e-6,
        direction_tolerance=1e-6,
        polarization_normal_tolerance=1e-6,
        polarization_fraction_tolerance=1e-6,
    ):
        self.wavelength_tolerance = wavelength_tolerance
        self.direction_tolerance = direction_tolerance
        self.polarization_normal_tolerance = polarization_normal_tolerance
        self.polarization_fraction_tolerance = polarization_fraction_tolerance

    def __call__(self, a, b):
        aw = a.get_wavelength()
        bw = b.get_wavelength()
        ad = matrix.col(a.get_sample_to_source_direction())
        bd = matrix.col(b.get_sample_to_source_direction())
        an = matrix.col(a.get_polarization_normal())
        bn = matrix.col(b.get_polarization_normal())
        af = a.get_polarization_fraction()
        bf = b.get_polarization_fraction()
        text = []
        if abs(aw - bw) > self.wavelength_tolerance:
            text.append(f" Wavelength: {aw:f}, {bw:f}")
        if abs(ad.angle(bd)) > self.direction_tolerance:
            text.append(" Direction: {}, {}".format(tuple(ad), tuple(bd)))
        if abs(an.angle(bn)) > self.polarization_normal_tolerance:
            text.append(" Polarization Normal: {}, {}".format(tuple(an), tuple(bn)))
        if abs(af - bf) > self.polarization_fraction_tolerance:
            text.append(f" Polarization Fraction: {af}, {bf}")
        if len(text) > 0:
            text = ["Beam:"] + text
        return text


class DetectorDiff:
    """A class to provide simple detector comparison"""

    def __init__(
        self, fast_axis_tolerance=1e-6, slow_axis_tolerance=1e-6, origin_tolerance=1e-6
    ):

        self.fast_axis_tolerance = fast_axis_tolerance
        self.slow_axis_tolerance = slow_axis_tolerance
        self.origin_tolerance = origin_tolerance

    def __call__(self, a, b):
        text = []
        if len(a) != len(b):
            text.append("Num Panels: %d, %d" % (len(a), len(b)))
        for i, (aa, bb) in enumerate(zip(a, b)):
            a_image_size = aa.get_image_size()
            b_image_size = bb.get_image_size()
            a_pixel_size = aa.get_pixel_size()
            b_pixel_size = bb.get_pixel_size()
            a_trusted_range = aa.get_trusted_range()
            b_trusted_range = bb.get_trusted_range()
            a_fast = aa.get_fast_axis()
            b_fast = bb.get_fast_axis()
            a_slow = aa.get_slow_axis()
            b_slow = bb.get_slow_axis()
            a_origin = aa.get_origin()
            b_origin = bb.get_origin()
            temp_text = []
            if not _all_equal(a_image_size, b_image_size):
                temp_text.append(f"  Image size: {a_image_size}, {b_image_size}")
            if not _all_approx_equal(a_pixel_size, b_pixel_size, 1e-7):
                temp_text.append(f"  Pixel size: {a_pixel_size}, {b_pixel_size}")
            if not _all_approx_equal(a_trusted_range, b_trusted_range, 1e-7):
                temp_text.append(
                    f"  Trusted Range: {a_trusted_range}, {b_trusted_range}"
                )
            if not _all_approx_equal(a_fast, b_fast, self.fast_axis_tolerance):
                temp_text.append(f"  Fast axis: {a_fast}, {b_fast}")
            if not _all_approx_equal(a_slow, b_slow, self.slow_axis_tolerance):
                temp_text.append(f"  Slow axis: {a_slow}, {b_slow}")
            if not _all_approx_equal(a_origin, b_origin, self.origin_tolerance):
                temp_text.append(f"  Origin: {a_origin}, {b_origin}")
            if len(temp_text) > 0:
                text.append(" panel %d:" % i)
                text.extend(temp_text)
        if len(text) > 0:
            text = ["Detector:"] + text
        return text


class GoniometerDiff:
    """A class to provide simple goniometer comparison"""

    def __init__(
        self,
        rotation_axis_tolerance=1e-6,
        fixed_rotation_tolerance=1e-6,
        setting_rotation_tolerance=1e-6,
    ):
        self.rotation_axis_tolerance = rotation_axis_tolerance
        self.fixed_rotation_tolerance = fixed_rotation_tolerance
        self.setting_rotation_tolerance = setting_rotation_tolerance

    def __call__(self, a, b):
        a_axis = matrix.col(a.get_rotation_axis())
        b_axis = matrix.col(b.get_rotation_axis())
        a_fixed = a.get_fixed_rotation()
        b_fixed = b.get_fixed_rotation()
        a_setting = a.get_setting_rotation()
        b_setting = b.get_setting_rotation()
        text = []
        if abs(a_axis.angle(b_axis)) > self.rotation_axis_tolerance:
            text.append(" Rotation axis: {}, {}".format(tuple(a_axis), tuple(b_axis)))
        if not _all_approx_equal(a_fixed, b_fixed, self.fixed_rotation_tolerance):
            text.append(f" Fixed rotation: {a_fixed}, {b_fixed}")
        if not _all_approx_equal(a_setting, b_setting, self.setting_rotation_tolerance):
            text.append(f" Setting rotation: {a_setting}, {b_setting}")
        if len(text) > 0:
            text = ["Goniometer:"] + text
        return text


class ScanDiff:
    """A class to provide scan comparison"""

    def __init__(self, scan_tolerance=1e-6):
        self.scan_tolerance = scan_tolerance

    def __call__(self, a, b):
        eps = self.scan_tolerance * abs(a.get_oscillation()[1])
        a_image_range = a.get_image_range()
        b_image_range = b.get_image_range()
        a_oscillation = a.get_oscillation()
        b_oscillation = b.get_oscillation()
        a_osc_range = a.get_oscillation_range()
        b_osc_range = b.get_oscillation_range()

        def mod_2pi(x):
            return x - 2 * math.pi * math.floor(x / (2 * math.pi))

        diff_2pi = abs(mod_2pi(a_osc_range[1]) - mod_2pi(b_osc_range[0]))
        diff_abs = abs(a_osc_range[1] - b_osc_range[0])
        text = []
        if not (a_image_range[1] + 1 == b_image_range[0]):
            text.append(f" Incompatible image range: {a_image_range}, {b_image_range}")
        if abs(a_oscillation[1] - b_oscillation[1]) > eps:
            text.append(f" Incompatible Oscillation: {a_oscillation}, {b_oscillation}")
        if min(diff_2pi, diff_abs) > eps * a.get_num_images():
            text.append(
                f" Incompatible Oscillation Range: {a_osc_range}, {b_osc_range}"
            )
        if len(text) > 0:
            text = ["Scan:"] + text
        return text


class SequenceDiff:
    def __init__(self, tolerance):

        if tolerance is None:
            self.b_diff = BeamDiff()
            self.d_diff = DetectorDiff()
            self.g_diff = GoniometerDiff()
            self.s_diff = ScanDiff()
        else:
            self.b_diff = BeamDiff(
                wavelength_tolerance=tolerance.beam.wavelength,
                direction_tolerance=tolerance.beam.direction,
                polarization_normal_tolerance=tolerance.beam.polarization_normal,
                polarization_fraction_tolerance=tolerance.beam.polarization_fraction,
            )

            self.d_diff = DetectorDiff(
                fast_axis_tolerance=tolerance.detector.fast_axis,
                slow_axis_tolerance=tolerance.detector.slow_axis,
                origin_tolerance=tolerance.detector.origin,
            )

            self.g_diff = GoniometerDiff(
                rotation_axis_tolerance=tolerance.goniometer.rotation_axis,
                fixed_rotation_tolerance=tolerance.goniometer.fixed_rotation,
                setting_rotation_tolerance=tolerance.goniometer.setting_rotation,
            )

            self.s_diff = ScanDiff(scan_tolerance=tolerance.scan.oscillation)

    def __call__(self, sequence1, sequence2):
        text = []
        text.extend(self.b_diff(sequence1.get_beam(), sequence2.get_beam()))
        text.extend(self.d_diff(sequence1.get_detector(), sequence2.get_detector()))
        text.extend(self.g_diff(sequence1.get_goniometer(), sequence2.get_goniometer()))
        text.extend(self.s_diff(sequence1.get_scan(), sequence2.get_scan()))
        return text
