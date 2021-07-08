import collections
import json
import logging
import os.path
import pickle
import warnings

import libtbx

import dxtbx.imageset
import dxtbx.model
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.format.image import ImageBool, ImageDouble
from dxtbx.format.Registry import get_format_class_for_file
from dxtbx.sequence_filenames import (
    locate_files_matching_template_string,
    template_string_number_index,
)
from dxtbx.serialize import load
from dxtbx.serialize.filename import resolve_path

logger = logging.getLogger(__name__)


class DataBlock:
    """High level container for blocks of sequences and imagesets."""

    def __init__(self, imagesets=None):
        """Instantiate from a list of imagesets."""

        warnings.warn(
            "Datablocks are deprecated; please use ExperimentLists instead",
            DeprecationWarning,
            stacklevel=2,
        )

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
            if isinstance(iset, dxtbx.imageset.RotImageSequence):
                yield iset

    def iter_stills(self):
        """Iterate over still groups."""
        for iset in self._imagesets:
            if not isinstance(iset, dxtbx.imageset.RotImageSequence):
                yield iset

    def _find_unique_items(self, item_name, filter_none=False):
        """Return a list of unique beams, detectors, ... in order.
        Optionally filter out None values (unless they came in via
        an RotImageSequence)."""
        items = {}
        for imageset in self._imagesets:
            getter_function = getattr(imageset, "get_" + item_name)
            if isinstance(imageset, dxtbx.imageset.RotImageSequence):
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
        result = {
            "__id__": "DataBlock",
            "imageset": [],
        }

        # Loop through all the imagesets
        for iset in self._imagesets:
            if isinstance(iset, dxtbx.imageset.RotImageSequence):
                if iset.reader().is_single_file_reader():
                    result["imageset"].append(
                        dict(
                            [
                                ("__id__", "RotImageSequence"),
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
                        dict(
                            [
                                ("__id__", "RotImageSequence"),
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
                imageset = {}
                if isinstance(iset, dxtbx.imageset.ImageGrid):
                    imageset["__id__"] = "ImageGrid"
                    imageset["grid_size"] = iset.get_grid_size()
                else:
                    imageset["__id__"] = "ImageSet"
                image_list = []
                for i in range(len(iset)):
                    image_dict = {
                        "filename": os.path.abspath(iset.get_path(i)),
                        "gain": abspath_or_none(iset.external_lookup.gain.filename),
                        "pedestal": abspath_or_none(
                            iset.external_lookup.pedestal.filename
                        ),
                        "dx": abspath_or_none(iset.external_lookup.dx.filename),
                        "dy": abspath_or_none(iset.external_lookup.dy.filename),
                    }
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


class InvalidDataBlockError(RuntimeError):
    """
    Indicates an error whilst validating the experiment list.

    This means that there is some structural problem that prevents the given data
    from representing a well-formed experiment list. This doesn't indicate e.g.
    some problem with the data or model consistency.
    """


def datablocks_from_dict(obj, check_format=True, directory=None):
    """Get the datablocks from the dictionary."""

    # If we have a list, extract for each dictionary in the list
    if isinstance(obj, list):
        return [datablocks_from_dict(dd, check_format, directory) for dd in obj]
    elif not isinstance(obj, dict):
        raise InvalidDataBlockError(
            "Unexpected datablock type {} instead of dict".format(type(obj))
        )
    # Make sure the id signature is correct
    if not obj.get("__id__") == "DataBlock":
        raise InvalidDataBlockError(
            "Expected __id__ 'DataBlock', but found {}".format(repr(obj.get("__id__")))
        )

    # Get the list of models
    blist = obj.get("beam", [])
    dlist = obj.get("detector", [])
    glist = obj.get("goniometer", [])
    slist = obj.get("scan", [])

    def load_models(obj):
        try:
            beam = dxtbx.model.BeamFactory.monochromatic_from_dict(blist[obj["beam"]])
        except Exception:
            beam = None
        try:
            dobj = dlist[obj["detector"]]
            detector = dxtbx.model.DetectorFactory.from_dict(dobj)
        except Exception:
            detector = None
        try:
            gonio = dxtbx.model.GoniometerFactory.from_dict(glist[obj["goniometer"]])
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
        if ident == "RotImageSequence" or ident == "ImageSweep":
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
                    imageset["dx"] = resolve_path(imageset["dx"], directory=directory)
                    iset.external_lookup.dx.filename = imageset["dx"]
                    with open(imageset["dx"], "rb") as infile:
                        iset.external_lookup.dx.data = ImageDouble(
                            pickle.load(infile, encoding="bytes")
                        )
                if "dy" in imageset and imageset["dy"] is not None:
                    imageset["dy"] = resolve_path(imageset["dy"], directory=directory)
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
            raise RuntimeError("expected ImageSet/RotImageSequence, got %s" % ident)

    return DataBlock(imagesets)


def datablocks_from_imagesets(imagesets):
    """Load a list of datablocks from imagesets."""

    datablocks = []
    if not isinstance(imagesets, list):
        imagesets = [imagesets]
    for imageset in imagesets:
        try:
            datablocks[-1].append(imageset)
        except (IndexError, AssertionError):
            datablocks.append(DataBlock([imageset]))
    return datablocks


class DataBlockFactory:
    """Class for creating DataBlock instances"""

    @staticmethod
    def from_args(
        args,
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
        unhandled=None,
        compare_beam=None,
        compare_detector=None,
        compare_goniometer=None,
        scan_tolerance=None,
        format_kwargs=None,
    ):
        """Create a list of data blocks from a list of directory or file names."""
        from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory

        expts = ExperimentListFactory.from_filenames(
            filenames,
            unhandled=unhandled,
            compare_beam=compare_beam,
            compare_detector=compare_detector,
            compare_goniometer=compare_goniometer,
            scan_tolerance=scan_tolerance,
            format_kwargs=format_kwargs,
        )
        datablocks = []

        # Datablocks can only contain imageset with a shared format class, therefore
        # group the experiments by format class
        format_groups = collections.defaultdict(ExperimentList)
        for expt in expts:
            format_groups[expt.imageset.get_format_class()].append(expt)

        # Generate a datablock for each format group
        for format_group, format_expts in format_groups.items():
            datablocks.extend(format_expts.to_datablocks())

        return datablocks

    @staticmethod
    def from_dict(obj, check_format=True, directory=None):
        """Create a datablock from a dictionary."""
        return datablocks_from_dict(obj, check_format, directory)

    @staticmethod
    def from_json(string, check_format=True, directory=None):
        """Decode a datablock from JSON string."""
        return DataBlockFactory.from_dict(
            json.loads(string),
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
        return datablocks_from_imagesets(imagesets)

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
