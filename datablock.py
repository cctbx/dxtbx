from __future__ import absolute_import, division, print_function

from builtins import range
import collections
import itertools
import json
import logging
import operator
import math
import os.path
from os.path import abspath, dirname, isdir, isfile, join, normpath, splitext
import warnings

import six.moves.cPickle as pickle

import dxtbx.imageset
import libtbx
import six
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.format.image import ImageBool, ImageDouble
from dxtbx.format.Registry import get_format_class_for_file
from dxtbx.model import BeamFactory, DetectorFactory, GoniometerFactory, ScanFactory
from dxtbx.serialize import load
from dxtbx.serialize.filename import resolve_path
from dxtbx.sweep_filenames import (
    locate_files_matching_template_string,
    template_regex,
    template_string_number_index,
)
from libtbx.utils import Sorry
from scitbx import matrix


logger = logging.getLogger(__name__)


class DataBlock(object):
    """ High level container for blocks of sweeps and imagesets. """

    def __init__(self, imagesets):
        """ Instantiate from a list of imagesets. """
        # Try to get a format class
        self._format_class = None
        self._imagesets = []

        # Load imagesets
        if imagesets is not None:
            for iset in imagesets:
                self.append(iset)

    def append(self, imageset):
        """ Add an imageset to the block. """
        if self._format_class is None:
            self._format_class = imageset.get_format_class()
        else:
            assert self._format_class == imageset.get_format_class()
        self._imagesets.append(imageset)

    def extend(self, datablock):
        """ Add two datablocks. """
        for iset in datablock:
            self.append(iset)

    def format_class(self):
        """ Return the format class. """
        return self._format_class

    def extract_stills(self):
        """ Extract all the still imagesets """
        return list(self.iter_stills())

    def extract_sweeps(self):
        """ Extract all the sweeps from the block. """
        return list(self.iter_sweeps())

    def extract_imagesets(self):
        """ Extract all imagesets. """
        return list(self._imagesets)

    def num_images(self):
        """ Get the number of images. """
        return sum(len(iset) for iset in self._imagesets)

    def __len__(self):
        """ The number of image sets. """
        return len(self._imagesets)

    def __eq__(self, rhs):
        """ Check if two blocks are the same. """
        return (
            self._format_class == rhs._format_class
            and self._imagesets == rhs._imagesets
        )

    __hash__ = None  # datablock objects are mutable and therefore unhashable

    def __ne__(self, rhs):
        """ Check if two blocks are not equal. """
        return not self.__eq__(rhs)

    def __iter__(self):
        """ Iterate through the imagesets. """
        for iset in self._imagesets:
            yield iset

    def iter_sweeps(self):
        """ Iterate over sweep groups. """
        for iset in self._imagesets:
            if isinstance(iset, dxtbx.imageset.ImageSweep):
                yield iset

    def iter_stills(self):
        """ Iterate over still groups. """
        for iset in self._imagesets:
            if not isinstance(iset, dxtbx.imageset.ImageSweep):
                yield iset

    def _find_unique_items(self, item_name, filter_none=False):
        """Return a list of unique beams, detectors, ... in order.
        Optionally filter out None values (unless they came in via
        an ImageSweep)."""

        # Use the keys of an ordered dictionary to guaranteeing uniqueness
        # while keeping the order. Could rewrite to use OrderedSet or wait
        # for Python 3.7 and use a regular dictionary
        # (3.7+: all dictionaries are guaranteed to be ordered)
        items = collections.OrderedDict()
        for imageset in self._imagesets:
            getter_function = getattr(imageset, "get_" + item_name)
            if isinstance(imageset, dxtbx.imageset.ImageSweep):
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
        """ Convert the datablock to a dictionary """

        def abspath_or_none(filename):
            if filename is None or filename == "":
                return None
            return abspath(filename)

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
            if isinstance(iset, dxtbx.imageset.ImageSweep):
                if iset.reader().is_single_file_reader():
                    result["imageset"].append(
                        collections.OrderedDict(
                            [
                                ("__id__", "ImageSweep"),
                                ("master", abspath(iset.reader().master_path())),
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
                                ("__id__", "ImageSweep"),
                                ("template", abspath(iset.get_template())),
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
                    image_dict["filename"] = abspath(iset.get_path(i))
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

        # Return the data block as a dictionary
        return result


class FormatChecker(object):
    """A helper class to speed up identifying the correct image format by first
    trying the last format that was used."""

    def __init__(self, verbose=None):
        """ Set the format class to none. """
        self._format_class = None
        if verbose is not None:
            warnings.warn(
                "The verbose parameter is deprecated.", DeprecationWarning, stacklevel=2
            )
        self._verbose = verbose

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


class DataBlockTemplateImporter(object):
    """ A class to import a datablock from a template. """

    def __init__(self, templates, verbose=None, **kwargs):
        """ Import the datablocks from the given templates. """
        if verbose is not None:
            warnings.warn(
                "The verbose parameter is deprecated.", DeprecationWarning, stacklevel=2
            )
        assert len(templates) > 0

        self.datablocks = []

        # A function to append or create a new datablock
        def append_to_datablocks(iset):
            try:
                self.datablocks[-1].append(iset)
            except Exception:
                self.datablocks.append(DataBlock([iset]))
            logger.debug("Added imageset to datablock %d", len(self.datablocks) - 1)

        # For each template do an import
        for template in templates:
            template = normpath(template)
            paths = sorted(locate_files_matching_template_string(template))
            logger.debug("The following files matched the template string:")
            if len(paths) > 0:
                for p in paths:
                    logger.debug(" %s", p)
            else:
                logger.debug(" No files found")

            # Check if we've matched any filenames
            if len(paths) == 0:
                raise Sorry('Template "%s" does not match any files' % template)

            # Get the format from the first image
            fmt = FormatChecker().find_format(paths[0])
            if fmt is None:
                raise Sorry("Image file %s format is unknown" % paths[0])
            elif fmt.ignore():
                raise Sorry("Image file %s format will be ignored" % paths[0])
            else:
                imageset = self._create_imageset(fmt, template, paths, **kwargs)
                append_to_datablocks(imageset)

    def _create_imageset(self, format_class, template, filenames, **kwargs):
        """ Create a multi file sweep or imageset. """
        # Get the image range
        index = slice(*template_string_number_index(template))
        first = int(filenames[0][index])
        last = int(filenames[-1][index])

        # Check all images in range are present - if allowed
        if not kwargs.get("allow_incomplete_sweeps", False):
            all_numbers = {int(f[index]) for f in filenames}
            missing = set(range(first, last + 1)) - all_numbers
            if missing:
                raise Sorry(
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

        # Create the sweep
        imageset = dxtbx.imageset.ImageSetFactory.make_sweep(
            template,
            list(range(*image_range)),
            format_class,
            b,
            d,
            g,
            s,
            format_kwargs=kwargs.get("format_kwargs"),
        )

        # Return the imageset
        return imageset


class DataBlockFilenameImporter(object):
    """ A class to import a datablock from image files. """

    def __init__(
        self,
        filenames,
        verbose=None,
        compare_beam=None,
        compare_detector=None,
        compare_goniometer=None,
        scan_tolerance=None,
        format_kwargs=None,
    ):
        """ Import the datablocks from the given filenames. """
        if verbose is not None:
            warnings.warn(
                "The verbose parameter is deprecated.", DeprecationWarning, stacklevel=2
            )

        # Init the datablock list
        self.unhandled = []
        self.datablocks = []

        # A function to append or create a new datablock
        def append_to_datablocks(iset):
            if self.datablocks:
                try:
                    self.datablocks[-1].append(iset)
                except AssertionError:
                    self.datablocks.append(DataBlock([iset]))
            else:
                self.datablocks.append(DataBlock([iset]))
            logger.debug("Added imageset to datablock %d", len(self.datablocks) - 1)

        # Iterate through groups of files by format class
        find_format = FormatChecker()
        for fmt, group in find_format.iter_groups(filenames):
            if fmt is None or fmt.ignore():
                self.unhandled.extend(group)
            elif issubclass(fmt, FormatMultiImage):
                for filename in group:
                    imageset = self._create_single_file_imageset(
                        fmt, filename, format_kwargs=format_kwargs
                    )
                    append_to_datablocks(imageset)
                    logger.debug("Loaded file: %s", filename)
            else:
                records = self._extract_file_metadata(
                    fmt,
                    group,
                    compare_beam,
                    compare_detector,
                    compare_goniometer,
                    scan_tolerance,
                    format_kwargs=format_kwargs,
                )
                for _, items in itertools.groupby(records, lambda r: r.group):
                    imageset = self._create_multi_file_imageset(
                        fmt, list(items), format_kwargs=format_kwargs
                    )
                    append_to_datablocks(imageset)

    def _extract_file_metadata(
        self,
        format_class,
        filenames,
        compare_beam=None,
        compare_detector=None,
        compare_goniometer=None,
        scan_tolerance=None,
        format_kwargs=None,
    ):
        """ Extract the file meta data in order to sort them. """
        # If no comparison functions are set
        if compare_beam is None:
            compare_beam = operator.__eq__
        if compare_detector is None:
            compare_detector = operator.__eq__
        if compare_goniometer is None:
            compare_goniometer = operator.__eq__
        if format_kwargs is None:
            format_kwargs = {}

        class Record(object):
            def __init__(
                self,
                beam=None,
                detector=None,
                goniometer=None,
                scan=None,
                template=None,
                filename=None,
                index=None,
                group=None,
            ):
                self.beam = beam
                self.detector = detector
                self.goniometer = goniometer
                self.scan = scan
                self.template = template
                self.filename = filename
                self.index = index
                self.group = group

        # Loop through all the filenames
        records = []
        group = 0
        for filename in filenames:

            # Read the image
            fmt = format_class(filename, **format_kwargs)

            # Get the meta data from the format
            try:
                b = fmt.get_beam()
            except Exception:
                b = None
            try:
                d = fmt.get_detector()
            except Exception:
                d = None
            try:
                g = fmt.get_goniometer()
            except Exception:
                g = None
            try:
                s = fmt.get_scan()
            except Exception:
                s = None

            # Get the template and index if possible
            if s is not None and abs(s.get_oscillation()[1]) > 0.0:
                template, index = template_regex(filename)
            else:
                template, index = None, None

            # Check the last record if available
            if len(records) > 0:
                last = records[-1]
                same = [False, False, False]
                if compare_beam(last.beam, b):
                    b = last.beam
                    same[0] = True
                if compare_detector(last.detector, d):
                    d = last.detector
                    same[1] = True
                if compare_goniometer(last.goniometer, g):
                    g = last.goniometer
                    same[2] = True

                # If the last was not a sweep then if the current is a sweep or none of
                # the models are the same, increment. Otherwise if the current is not a
                # sweep or if both sweeps don't share all models increment. If both
                # share, models try scan and increment if exception.
                if last.template is None:
                    if template is not None or not any(same):
                        group += 1
                else:
                    if template is None or not all(same):
                        group += 1
                    else:
                        try:
                            if scan_tolerance is None:
                                last.scan.append(s)
                            else:
                                last.scan.append(s, scan_tolerance=scan_tolerance)
                            last.index = index
                            continue
                        except Exception:
                            group += 1

            # Add a record
            records.append(
                Record(
                    beam=b,
                    detector=d,
                    goniometer=g,
                    scan=s,
                    template=template,
                    filename=filename,
                    index=index,
                    group=group,
                )
            )

        # Return the records
        return records

    def _create_multi_file_imageset(self, format_class, records, format_kwargs=None):
        """ Create a multi file sweep or imageset. """

        # Make either an imageset or sweep
        if len(records) == 1 and records[0].template is not None:

            # Get the image range
            image_range = records[0].scan.get_image_range()
            image_range = (image_range[0], image_range[1] + 1)

            # Create the sweep
            imageset = dxtbx.imageset.ImageSetFactory.make_sweep(
                abspath(records[0].template),
                list(range(*image_range)),
                format_class,
                records[0].beam,
                records[0].detector,
                records[0].goniometer,
                records[0].scan,
                format_kwargs=format_kwargs,
            )

        else:

            # Get the filenames
            filenames = []
            for r in records:
                assert r.template is None
                filenames.append(r.filename)

            # make an imageset
            imageset = dxtbx.imageset.ImageSetFactory.make_imageset(
                list(map(abspath, filenames)), format_class, format_kwargs=format_kwargs
            )
            for i, r in enumerate(records):
                imageset.set_beam(r.beam, i)
                imageset.set_detector(r.detector, i)
                imageset.set_goniometer(r.goniometer, i)
                imageset.set_scan(r.scan, i)

        # Return the imageset
        return imageset

    def _create_single_file_imageset(self, format_class, filename, format_kwargs=None):
        """ Create an imageset from a multi image file. """
        if format_kwargs is None:
            format_kwargs = {}
        return format_class.get_imageset(abspath(filename), format_kwargs=format_kwargs)


class InvalidDataBlockError(RuntimeError):
    """
    Indicates an error whilst validating the experiment list.

    This means that there is some structural problem that prevents the given data
    from representing a well-formed experiment list. This doesn't indicate e.g.
    some problem with the data or model consistency.
    """


class DataBlockDictImporter(object):
    """ A class to import a datablock from dictionary. """

    def __init__(self, obj, check_format=True, directory=None):
        """ Get the datablocks from the dictionary. """
        self.datablocks = self._load_datablocks(obj, check_format, directory)

    def _load_datablocks(self, obj, check_format=True, directory=None):
        """ Create the datablock from a dictionary. """

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
                beam = BeamFactory.from_dict(blist[obj["beam"]])
            except Exception:
                beam = None
            try:
                dobj = dlist[obj["detector"]]
                detector = DetectorFactory.from_dict(dobj)
            except Exception:
                detector = None
            try:
                gonio = GoniometerFactory.from_dict(glist[obj["goniometer"]])
            except Exception:
                gonio = None
            try:
                scan = ScanFactory.from_dict(slist[obj["scan"]])
            except Exception:
                scan = None
            return beam, detector, gonio, scan

        if six.PY3:
            pickle_parameters = {"encoding": "bytes"}
        else:
            pickle_parameters = {}

        # Loop through all the imagesets
        imagesets = []
        for imageset in obj["imageset"]:
            ident = imageset["__id__"]
            if "params" in imageset:
                format_kwargs = imageset["params"]
            else:
                format_kwargs = {}
            if ident == "ImageSweep":
                beam, detector, gonio, scan = load_models(imageset)
                if "template" in imageset:
                    template = resolve_path(imageset["template"], directory=directory)
                    i0, i1 = scan.get_image_range()
                    iset = dxtbx.imageset.ImageSetFactory.make_sweep(
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
                                    pickle.load(infile, **pickle_parameters)
                                )
                    if "gain" in imageset and imageset["gain"] is not None:
                        imageset["gain"] = resolve_path(
                            imageset["gain"], directory=directory
                        )
                        iset.external_lookup.gain.filename = imageset["gain"]
                        if check_format:
                            with open(imageset["gain"], "rb") as infile:
                                iset.external_lookup.gain.data = ImageDouble(
                                    pickle.load(infile, **pickle_parameters)
                                )
                    if "pedestal" in imageset and imageset["pedestal"] is not None:
                        imageset["pedestal"] = resolve_path(
                            imageset["pedestal"], directory=directory
                        )
                        iset.external_lookup.pedestal.filename = imageset["pedestal"]
                        if check_format:
                            with open(imageset["pedestal"], "rb") as infile:
                                iset.external_lookup.pedestal.data = ImageDouble(
                                    pickle.load(infile, **pickle_parameters)
                                )
                    if "dx" in imageset and imageset["dx"] is not None:
                        imageset["dx"] = resolve_path(
                            imageset["dx"], directory=directory
                        )
                        iset.external_lookup.dx.filename = imageset["dx"]
                        with open(imageset["dx"], "rb") as infile:
                            iset.external_lookup.dx.data = ImageDouble(
                                pickle.load(infile, **pickle_parameters)
                            )
                    if "dy" in imageset and imageset["dy"] is not None:
                        imageset["dy"] = resolve_path(
                            imageset["dy"], directory=directory
                        )
                        iset.external_lookup.dy.filename = imageset["dy"]
                        with open(imageset["dy"], "rb") as infile:
                            iset.external_lookup.dy.data = ImageDouble(
                                pickle.load(infile, **pickle_parameters)
                            )
                    iset.update_detector_px_mm_data()
                elif "master" in imageset:
                    template = resolve_path(imageset["master"], directory=directory)
                    i0, i1 = scan.get_image_range()
                    if not check_format:
                        format_class = FormatMultiImage
                    else:
                        format_class = None
                    iset = dxtbx.imageset.ImageSetFactory.make_sweep(
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
                                    pickle.load(infile, **pickle_parameters)
                                )
                    if "gain" in imageset and imageset["gain"] is not None:
                        imageset["gain"] = resolve_path(imageset["gain"], directory)
                        iset.external_lookup.gain.filename = imageset["gain"]
                        if check_format:
                            with open(imageset["gain"], "rb") as infile:
                                iset.external_lookup.gain.data = ImageDouble(
                                    pickle.load(infile, **pickle_parameters)
                                )
                    if "pedestal" in imageset and imageset["pedestal"] is not None:
                        imageset["pedestal"] = resolve_path(
                            imageset["pedestal"], directory
                        )
                        iset.external_lookup.pedestal.filename = imageset["pedestal"]
                        if check_format:
                            with open(imageset["pedestal"], "rb") as infile:
                                iset.external_lookup.pedestal.data = ImageDouble(
                                    pickle.load(infile, **pickle_parameters)
                                )
                    if "dx" in imageset and imageset["dx"] is not None:
                        imageset["dx"] = resolve_path(imageset["dx"], directory)
                        iset.external_lookup.dx.filename = imageset["dx"]
                        with open(imageset["dx"], "rb") as infile:
                            iset.external_lookup.dx.data = ImageDouble(
                                pickle.load(infile, **pickle_parameters)
                            )
                    if "dy" in imageset and imageset["dy"] is not None:
                        imageset["dy"] = resolve_path(imageset["dy"], directory)
                        iset.external_lookup.dy.filename = imageset["dy"]
                        with open(imageset["dy"], "rb") as infile:
                            iset.external_lookup.dy.data = ImageDouble(
                                pickle.load(infile, **pickle_parameters)
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
                                pickle.load(infile, **pickle_parameters)
                            )
                if "gain" in imageset and imageset["gain"] is not None:
                    imageset["gain"] = resolve_path(imageset["gain"], directory)
                    iset.external_lookup.gain.filename = imageset["gain"]
                    if check_format:
                        with open(imageset["gain"], "rb") as infile:
                            iset.external_lookup.gain.data = ImageDouble(
                                pickle.load(infile, **pickle_parameters)
                            )
                if "pedestal" in imageset and imageset["pedestal"] is not None:
                    imageset["pedestal"] = resolve_path(imageset["pedestal"], directory)
                    iset.external_lookup.pedestal.filename = imageset["pedestal"]
                    if check_format:
                        with open(imageset["pedestal"], "rb") as infile:
                            iset.external_lookup.pedestal.data = ImageDouble(
                                pickle.load(infile, **pickle_parameters)
                            )
                if "dx" in imageset and imageset["dx"] is not None:
                    imageset["dx"] = resolve_path(imageset["dx"], directory)
                    iset.external_lookup.dx.filename = imageset["dx"]
                    with open(imageset["dx"], "rb") as infile:
                        iset.external_lookup.dx.data = ImageDouble(
                            pickle.load(infile, **pickle_parameters)
                        )
                if "dy" in imageset and imageset["dy"] is not None:
                    imageset["dy"] = resolve_path(imageset["dy"], directory)
                    iset.external_lookup.dy.filename = imageset["dy"]
                    with open(imageset["dy"], "rb") as infile:
                        iset.external_lookup.dy.data = ImageDouble(
                            pickle.load(infile, **pickle_parameters)
                        )
                    iset.update_detector_px_mm_data()
                imagesets.append(iset)
            else:
                raise RuntimeError("expected ImageSet/ImageSweep, got %s" % ident)

        # Return the datablock
        return DataBlock(imagesets)


class DataBlockImageSetImporter(object):
    """ A class to import a datablock from imagesets. """

    def __init__(self, imagesets):
        """ Load a list of datablocks from imagesets. """
        self.datablocks = []
        if not isinstance(imagesets, list):
            imagesets = [imagesets]
        for imageset in imagesets:
            try:
                self.datablocks[-1].append(imageset)
            except (IndexError, AssertionError):
                self.datablocks.append(DataBlock([imageset]))


class DataBlockFactory(object):
    """ Class for creating DataBlock instances"""

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
        """ Try to load datablocks from any recognized format. """

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

        # Return the datablocks
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
        """ Create a list of data blocks from a list of directory or file names. """
        filelist = []
        for f in sorted(filenames):
            if isfile(f):
                filelist.append(f)
            elif isdir(f):
                subdir = sorted(
                    join(f, sf) for sf in os.listdir(f) if isfile(join(f, sf))
                )
                filelist.extend(subdir)
                logger.debug("Added %d files from %s", len(subdir), f)
            else:
                logger.debug(
                    "Could not import %s: not a valid file or directory name", f
                )
                if unhandled is not None:
                    unhandled.append(f)

        importer = DataBlockFilenameImporter(
            filelist,
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
        """ Create a datablock from a dictionary. """
        importer = DataBlockDictImporter(obj, check_format, directory)
        return importer.datablocks

    @staticmethod
    def from_json(string, check_format=True, directory=None):
        """ Decode a datablock from JSON string. """
        from dxtbx.serialize.load import _decode_dict

        return DataBlockFactory.from_dict(
            json.loads(string, object_hook=_decode_dict),
            check_format=check_format,
            directory=directory,
        )

    @staticmethod
    def from_json_file(filename, check_format=True):
        """ Decode a datablock from a JSON file. """
        filename = abspath(filename)
        directory = dirname(filename)
        with open(filename, "r") as infile:
            return DataBlockFactory.from_json(
                infile.read(), check_format=check_format, directory=directory
            )

    @staticmethod
    def from_pickle_file(filename):
        """ Decode a datablock from a pickle file. """
        with open(filename, "rb") as infile:
            obj = pickle.load(infile)
            if isinstance(obj, list):
                assert all(isinstance(db, DataBlock) for db in obj)
            else:
                assert isinstance(obj, DataBlock)
            return obj

    @staticmethod
    def from_imageset(imagesets):
        """ Load a datablock from a list of imagesets. """
        importer = DataBlockImageSetImporter(imagesets)
        return importer.datablocks

    @staticmethod
    def from_imageset_json_file(filename):
        """ Load a datablock from a sweep file. """
        # Load the imageset and create a datablock from the filenames
        imageset = load.imageset(filename)
        return DataBlockFactory.from_imageset(imageset)

    @staticmethod
    def from_serialized_format(filename, check_format=True):
        """ Load a datablock from serialized formats. """

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
        """ Function to instantiate data block from in memory imageset. """
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


class DataBlockDumper(object):
    """ Class to help in dumping datablock objects. """

    def __init__(self, datablocks):
        """ Initialise the list of data blocks. """
        if isinstance(datablocks, DataBlock):
            self._datablocks = [datablocks]
        else:
            self._datablocks = datablocks

    def as_json(self, filename=None, compact=False):
        """ Dump datablock as json. """

        dictionary = [db.to_dict() for db in self._datablocks]
        if compact:
            json.dump(
                dictionary,
                open(filename, "w"),
                separators=(",", ":"),
                ensure_ascii=True,
                cls=AutoEncoder,
            )
        else:
            json.dump(
                dictionary,
                open(filename, "w"),
                indent=2,
                ensure_ascii=True,
                cls=AutoEncoder,
            )

    def as_pickle(self, filename=None, **kwargs):
        """ Dump datablock as pickle. """

        # Get the pickle string
        text = pickle.dumps(self._datablocks, protocol=pickle.HIGHEST_PROTOCOL)

        # Write the file
        if filename is not None:
            with open(filename, "wb") as outfile:
                outfile.write(text)
        else:
            return text

    def as_file(self, filename, **kwargs):
        """ Dump datablocks as file. """
        ext = splitext(filename)[1]
        j_ext = [".json"]
        p_ext = [".p", ".pkl", ".pickle"]
        if ext.lower() in j_ext:
            return self.as_json(filename, **kwargs)
        elif ext.lower() in p_ext:
            return self.as_pickle(filename, **kwargs)
        else:
            ext_str = "|".join(j_ext + p_ext)
            raise RuntimeError("expected extension {%s}, got %s" % (ext_str, ext))


class BeamComparison(object):
    """
    A class to provide simple beam comparison

    """

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


class DetectorComparison(object):
    """
    A class to provide simple detector comparison

    """

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


class GoniometerComparison(object):
    """
    A class to provide simple goniometer comparison

    """

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


class BeamDiff(object):
    """
    A class to provide simple beam comparison

    """

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
        ad = matrix.col(a.get_direction())
        bd = matrix.col(b.get_direction())
        an = matrix.col(a.get_polarization_normal())
        bn = matrix.col(b.get_polarization_normal())
        af = a.get_polarization_fraction()
        bf = b.get_polarization_fraction()
        text = []
        if abs(aw - bw) > self.wavelength_tolerance:
            text.append(" Wavelength: %f, %f" % (aw, bw))
        if abs(ad.angle(bd)) > self.direction_tolerance:
            text.append(" Direction: %s, %s" % (tuple(ad), tuple(bd)))
        if abs(an.angle(bn)) > self.polarization_normal_tolerance:
            text.append(" Polarization Normal: %s, %s" % (tuple(an), tuple(bn)))
        if abs(af - bf) > self.polarization_fraction_tolerance:
            text.append(" Polarization Fraction: %s, %s" % (af, bf))
        if len(text) > 0:
            text = ["Beam:"] + text
        return text


class DetectorDiff(object):
    """
    A class to provide simple detector comparison

    """

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
                temp_text.append("  Image size: %s, %s" % (a_image_size, b_image_size))
            if not _all_approx_equal(a_pixel_size, b_pixel_size, 1e-7):
                temp_text.append("  Pixel size: %s, %s" % (a_pixel_size, b_pixel_size))
            if not _all_approx_equal(a_trusted_range, b_trusted_range, 1e-7):
                temp_text.append(
                    "  Trusted Range: %s, %s" % (a_trusted_range, b_trusted_range)
                )
            if not _all_approx_equal(a_fast, b_fast, self.fast_axis_tolerance):
                temp_text.append("  Fast axis: %s, %s" % (a_fast, b_fast))
            if not _all_approx_equal(a_slow, b_slow, self.slow_axis_tolerance):
                temp_text.append("  Slow axis: %s, %s" % (a_slow, b_slow))
            if not _all_approx_equal(a_origin, b_origin, self.origin_tolerance):
                temp_text.append("  Origin: %s, %s" % (a_origin, b_origin))
            if len(temp_text) > 0:
                text.append(" panel %d:" % i)
                text.extend(temp_text)
        if len(text) > 0:
            text = ["Detector:"] + text
        return text


class GoniometerDiff(object):
    """
    A class to provide simple goniometer comparison

    """

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
            text.append(" Rotation axis: %s, %s" % (tuple(a_axis), tuple(b_axis)))
        if not _all_approx_equal(a_fixed, b_fixed, self.fixed_rotation_tolerance):
            text.append(" Fixed rotation: %s, %s" % (a_fixed, b_fixed))
        if not _all_approx_equal(a_setting, b_setting, self.setting_rotation_tolerance):
            text.append(" Setting rotation: %s, %s" % (a_setting, b_setting))
        if len(text) > 0:
            text = ["Goniometer:"] + text
        return text


class ScanDiff(object):
    """
    A class to provide scan comparison

    """

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
            text.append(
                " Incompatible image range: %s, %s" % (a_image_range, b_image_range)
            )
        if abs(a_oscillation[1] - b_oscillation[1]) > eps:
            text.append(
                " Incompatible Oscillation: %s, %s" % (a_oscillation, b_oscillation)
            )
        if min(diff_2pi, diff_abs) > eps * a.get_num_images():
            text.append(
                " Incompatible Oscillation Range: %s, %s" % (a_osc_range, b_osc_range)
            )
        if len(text) > 0:
            text = ["Scan:"] + text
        return text


class SweepDiff(object):
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

    def __call__(self, sweep1, sweep2):
        text = []
        text.extend(self.b_diff(sweep1.get_beam(), sweep2.get_beam()))
        text.extend(self.d_diff(sweep1.get_detector(), sweep2.get_detector()))
        text.extend(self.g_diff(sweep1.get_goniometer(), sweep2.get_goniometer()))
        text.extend(self.s_diff(sweep1.get_scan(), sweep2.get_scan()))
        return text
