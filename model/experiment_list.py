import collections
import copy
import errno
import itertools
import json
import logging
import operator
import os
import pickle
import warnings

import pkg_resources

import dxtbx.datablock
from dxtbx.datablock import (
    BeamComparison,
    DataBlockFactory,
    DataBlockTemplateImporter,
    DetectorComparison,
    GoniometerComparison,
)
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.format.image import ImageBool, ImageDouble
from dxtbx.imageset import ImageGrid, ImageSequence, ImageSet, ImageSetFactory
from dxtbx.model import (
    BeamFactory,
    CrystalFactory,
    DetectorFactory,
    Experiment,
    ExperimentList,
    GoniometerFactory,
    ProfileModelFactory,
    ScanFactory,
)
from dxtbx.sequence_filenames import template_image_range, template_regex
from dxtbx.serialize import xds
from dxtbx.serialize.filename import resolve_path
from dxtbx.util import get_url_scheme

try:
    from typing import (
        Any,
        Callable,
        Dict,
        Generator,
        Iterable,
        List,
        Optional,
        Tuple,
        Type,
    )
except ImportError:
    pass

__all__ = [
    "BeamComparison",
    "DetectorComparison",
    "ExperimentListFactory",
    "GoniometerComparison",
]


logger = logging.getLogger(__name__)


class InvalidExperimentListError(RuntimeError):
    """
    Indicates an error whilst validating the experiment list.

    This means that there is some structural problem that prevents the given data
    from representing a well-formed experiment list. This doesn't indicate e.g.
    some problem with the data or model consistency.
    """


class ExperimentListDict(object):
    """A helper class for serializing the experiment list to dictionary (needed
    to save the experiment list to JSON format."""

    def __init__(self, obj, check_format=True, directory=None):
        """Initialise. Copy the dictionary."""
        # Basic check: This is a dict-like object. This can happen if e.g. we
        # were passed a DataBlock list instead of an ExperimentList dictionary
        if isinstance(obj, list) or not hasattr(obj, "get"):
            raise InvalidExperimentListError(
                "Expected dictionary, not {}".format(type(obj))
            )

        self._obj = copy.deepcopy(obj)
        self._check_format = check_format
        self._directory = directory

        # If this doesn't claim to be an ExperimentList, don't even try
        if self._obj.get("__id__") != "ExperimentList":
            raise InvalidExperimentListError(
                "Expected __id__ 'ExperimentList', but found {}".format(
                    repr(self._obj.get("__id__"))
                )
            )

        # Extract lists of models referenced by experiments
        # Go through all the imagesets and make sure the dictionary
        # references by an index rather than a file path.
        self._lookups = {
            model: self._extract_models(model, function)
            for model, function in (
                ("beam", BeamFactory.from_dict),
                ("detector", DetectorFactory.from_dict),
                ("goniometer", GoniometerFactory.from_dict),
                ("scan", ScanFactory.from_dict),
                ("crystal", CrystalFactory.from_dict),
                ("profile", ProfileModelFactory.from_dict),
                ("imageset", lambda x: x),
                ("scaling_model", self._scaling_model_from_dict),
            )
        }

    def _extract_models(self, name, from_dict):
        """
        Helper function. Extract the models.

        if name == imageset: Extract imageset objects from the source.

        This function does resolving of an (old) method of imageset lookup
        e.g. it was valid to have a string as the imageset value in an
        experiment instead of an int - in which case the imageset was
        loaded from the named file in the target directory.

        If any experiments point to a file in this way, the imageset is
        loaded and the experiment is rewritted with an integer pointing
        to the new ImageSet in the returned list.

        Returns:
                The ordered list of serialized-ImageSet dictionaries
                that the Experiment list points to.
        """

        # Extract all the model list
        mlist = self._obj.get(name, [])

        # Convert the model from dictionary to concreate
        # python class for the model.
        mlist = [from_dict(d) for d in mlist]

        # Dictionaries for file mappings
        mmap = {}

        # For each experiment, check the model is not specified by
        # a path, if it is then get the dictionary of the model
        # and insert it into the list. Replace the path reference
        # with an index
        for eobj in self._obj["experiment"]:
            value = eobj.get(name)
            if value is None:
                continue
            elif isinstance(value, str):
                if value not in mmap:
                    mmap[value] = len(mlist)
                    mlist.append(
                        from_dict(_experimentlist_from_file(value, self._directory))
                    )
                eobj[name] = mmap[value]
            elif not isinstance(value, int):
                raise TypeError("expected int or str, got %s" % type(value))

        return mlist

    def _load_pickle_path(self, imageset_data, param):
        # type: (Dict, str) -> Tuple[Optional[str], Any]
        """
        Read a filename from an imageset dict and load if required.

        Args:
            imageset_data: The dictionary holding imageset information
            param: The key name to lookup in the imageset dictionary

        Returns:
            A tuple of (filename, data) where data has been loaded from
            the pickle file. If there is no key entry then (None, None)
            is returned. If the configuration parameter check_format is
            False then (filename, None) will be returned.
        """
        if param not in imageset_data:
            return "", None

        filename = resolve_path(imageset_data[param], directory=self._directory)
        if self._check_format and filename:
            with open(filename, "rb") as fh:
                return filename, pickle.load(fh, encoding="bytes")

        return filename or "", None

    def _imageset_from_imageset_data(self, imageset_data, models):
        """Make an imageset from imageset_data - help with refactor decode."""
        assert imageset_data is not None
        if "params" in imageset_data:
            format_kwargs = imageset_data["params"]
        else:
            format_kwargs = {}

        beam = models["beam"]
        detector = models["detector"]
        goniometer = models["goniometer"]
        scan = models["scan"]

        # Load the external lookup data
        mask_filename, mask = self._load_pickle_path(imageset_data, "mask")
        gain_filename, gain = self._load_pickle_path(imageset_data, "gain")
        pedestal_filename, pedestal = self._load_pickle_path(imageset_data, "pedestal")
        dx_filename, dx = self._load_pickle_path(imageset_data, "dx")
        dy_filename, dy = self._load_pickle_path(imageset_data, "dy")

        if imageset_data["__id__"] == "ImageSet":
            imageset = self._make_stills(imageset_data, format_kwargs=format_kwargs)
        elif imageset_data["__id__"] == "ImageGrid":
            imageset = self._make_grid(imageset_data, format_kwargs=format_kwargs)
        elif (
            imageset_data["__id__"] == "ImageSequence"
            or imageset_data["__id__"] == "ImageSweep"
        ):
            imageset = self._make_sequence(
                imageset_data,
                beam=beam,
                detector=detector,
                goniometer=goniometer,
                scan=scan,
                format_kwargs=format_kwargs,
            )
        elif imageset_data["__id__"] == "MemImageSet":
            imageset = self._make_mem_imageset(imageset_data)
        else:
            raise RuntimeError("Unknown imageset type")

        if imageset is not None:
            # Set the external lookup
            if mask is None:
                mask = ImageBool()
            else:
                mask = ImageBool(mask)
            if gain is None:
                gain = ImageDouble()
            else:
                gain = ImageDouble(gain)
            if pedestal is None:
                pedestal = ImageDouble()
            else:
                pedestal = ImageDouble(pedestal)
            if dx is None:
                dx = ImageDouble()
            else:
                dx = ImageDouble(dx)
            if dy is None:
                dy = ImageDouble()
            else:
                dy = ImageDouble(dy)

            if not imageset.external_lookup.mask.data.empty():
                if not mask.empty():
                    mask = tuple(m.data() for m in mask)
                    for m1, m2 in zip(mask, imageset.external_lookup.mask.data):
                        m1 &= m2.data()
                    imageset.external_lookup.mask.data = ImageBool(mask)
            else:
                imageset.external_lookup.mask.data = mask
            imageset.external_lookup.mask.filename = mask_filename
            imageset.external_lookup.gain.data = gain
            imageset.external_lookup.gain.filename = gain_filename
            imageset.external_lookup.pedestal.data = pedestal
            imageset.external_lookup.pedestal.filename = pedestal_filename
            imageset.external_lookup.dx.data = dx
            imageset.external_lookup.dx.filename = dx_filename
            imageset.external_lookup.dy.data = dy
            imageset.external_lookup.dy.filename = dy_filename

            # Update the imageset models
            if isinstance(imageset, ImageSequence):
                imageset.set_beam(beam)
                imageset.set_detector(detector)
                imageset.set_goniometer(goniometer)
                imageset.set_scan(scan)
            elif isinstance(imageset, (ImageSet, ImageGrid)):
                for i in range(len(imageset)):
                    imageset.set_beam(beam, i)
                    imageset.set_detector(detector, i)
                    imageset.set_goniometer(goniometer, i)
                    imageset.set_scan(scan, i)

            imageset.update_detector_px_mm_data()

        return imageset

    def decode(self):
        """Decode the dictionary into a list of experiments."""
        # Extract all the experiments - first find all scans belonging to
        # same imageset

        eobj_scan = {}

        for eobj in self._obj["experiment"]:
            if self._lookup_model("imageset", eobj) is None:
                continue
            imageset_ref = eobj.get("imageset")
            scan = self._lookup_model("scan", eobj)

            if imageset_ref in eobj_scan:
                # if there is no scan, or scan is identical, move on, else
                # make a scan which encompasses both scans
                if not scan or scan == eobj_scan[imageset_ref]:
                    continue
                i = eobj_scan[imageset_ref].get_image_range()
                j = scan.get_image_range()
                if i[1] + 1 == j[0]:
                    eobj_scan[imageset_ref] += scan
                else:
                    # make a new bigger scan
                    o = eobj_scan[imageset_ref].get_oscillation()
                    s = scan.get_oscillation()
                    assert o[1] == s[1]
                    scan = copy.deepcopy(scan)
                    scan.set_image_range((min(i[0], j[0]), max(i[1], j[1])))
                    scan.set_oscillation((min(o[0], s[0]), o[1]))
                    eobj_scan[imageset_ref] = scan
            else:
                eobj_scan[imageset_ref] = copy.deepcopy(scan)

        # Map of imageset/scan pairs
        imagesets = {}

        # For every experiment, use the given input to create
        # a sensible experiment.
        el = ExperimentList()
        for eobj in self._obj["experiment"]:

            # Get the models
            identifier = eobj.get("identifier", "")
            beam = self._lookup_model("beam", eobj)
            detector = self._lookup_model("detector", eobj)
            goniometer = self._lookup_model("goniometer", eobj)
            scan = self._lookup_model("scan", eobj)
            crystal = self._lookup_model("crystal", eobj)
            profile = self._lookup_model("profile", eobj)
            scaling_model = self._lookup_model("scaling_model", eobj)

            models = {
                "beam": beam,
                "detector": detector,
                "goniometer": goniometer,
                "scan": scan,
                "crystal": crystal,
                "profile": profile,
                "scaling_model": scaling_model,
            }

            imageset_ref = eobj.get("imageset")

            # If not already cached, load this imageset
            if imageset_ref not in imagesets:
                imageset_data = self._lookup_model("imageset", eobj)
                if imageset_data is not None:
                    # Create the imageset from the input data
                    models["scan"] = eobj_scan[imageset_ref]
                    imageset = self._imageset_from_imageset_data(imageset_data, models)
                    imagesets[imageset_ref] = imageset
                else:
                    # Even if we have an empty entry, this counts as a load
                    imagesets[imageset_ref] = None

            # Append the experiment
            el.append(
                Experiment(
                    imageset=imagesets[imageset_ref],
                    beam=beam,
                    detector=detector,
                    goniometer=goniometer,
                    scan=scan,
                    crystal=crystal,
                    profile=profile,
                    scaling_model=scaling_model,
                    identifier=identifier,
                )
            )

        return el

    def _make_mem_imageset(self, imageset):
        """Can't make a mem imageset from dict."""
        return None

    def _make_stills(self, imageset, format_kwargs=None):
        """Make a still imageset."""
        filenames = [
            resolve_path(p, directory=self._directory) if not get_url_scheme(p) else p
            for p in imageset["images"]
        ]
        indices = None
        if "single_file_indices" in imageset:
            indices = imageset["single_file_indices"]
            assert len(indices) == len(filenames)
        return ImageSetFactory.make_imageset(
            filenames,
            None,
            check_format=self._check_format,
            single_file_indices=indices,
            format_kwargs=format_kwargs,
        )

    def _make_grid(self, imageset, format_kwargs=None):
        """Make a still imageset."""
        grid_size = imageset["grid_size"]
        return ImageGrid.from_imageset(
            self._make_stills(imageset, format_kwargs=format_kwargs), grid_size
        )

    def _make_sequence(
        self,
        imageset,
        beam=None,
        detector=None,
        goniometer=None,
        scan=None,
        format_kwargs=None,
    ):
        """Make an image sequence."""
        # Get the template format
        template = resolve_path(imageset["template"], directory=self._directory)

        # Get the number of images (if no scan is given we'll try
        # to find all the images matching the template
        if scan is None:
            i0, i1 = template_image_range(template)
        else:
            i0, i1 = scan.get_image_range()

        format_class = None
        if self._check_format is False:
            if "single_file_indices" in imageset:
                format_class = FormatMultiImage

        # Make a sequence from the input data
        return ImageSetFactory.make_sequence(
            template,
            list(range(i0, i1 + 1)),
            format_class=format_class,
            check_format=self._check_format,
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            format_kwargs=format_kwargs,
        )

    def _lookup_model(self, name, experiment_dict):
        """
        Find a model by looking up its index from a dictionary

        Args:
            name (str): The model name e.g. 'beam', 'detector'
            experiment_dict (Dict[str, int]):
                The experiment dictionary. experiment_dict[name] must
                exist and be not None to retrieve a model. If this key
                exists, then there *must* be an item with this index
                in the ExperimentListDict internal model stores.

        Returns:
            Optional[Any]:
                A model by looking up the index pointed to by
                experiment_dict[name]. If not present or empty,
                then None is returned.
        """
        if experiment_dict.get(name) is None:
            return None
        return self._lookups[name][experiment_dict[name]]

    @staticmethod
    def _scaling_model_from_dict(obj):
        """Get the scaling model from a dictionary."""
        for entry_point in pkg_resources.iter_entry_points("dxtbx.scaling_model_ext"):
            if entry_point.name == obj["__id__"]:
                return entry_point.load().from_dict(obj)


def _experimentlist_from_file(filename, directory=None):
    """Load a model dictionary from a file."""
    filename = resolve_path(filename, directory=directory)
    try:
        with open(filename, "r") as infile:
            return json.load(infile)
    except IOError:
        raise IOError("unable to read file, %s" % filename)


class ExperimentListFactory(object):
    """A class to help instantiate experiment lists."""

    @staticmethod
    def from_args(args, unhandled=None):
        """Try to load serialised experiments from any recognised format."""

        # Create a list for unhandled arguments
        if unhandled is None:
            unhandled = []

        experiments = ExperimentList()

        # Try to load from serialized formats
        for filename in args:
            try:
                experiments.extend(
                    ExperimentListFactory.from_serialized_format(filename)
                )
                logger.debug(f"Loaded experiments from {filename}")
            except Exception as e:
                logger.debug(f"Could not load experiments from {filename}: {e}")
                unhandled.append(filename)

        return experiments

    @staticmethod
    def from_filenames(
        filenames,
        unhandled=None,
        compare_beam=None,
        compare_detector=None,
        compare_goniometer=None,
        scan_tolerance=None,
        format_kwargs=None,
        load_models=True,
    ):
        """Create a list of data blocks from a list of directory or file names."""
        experiments = ExperimentList()

        # Process each file given by this path list
        to_process = _openingpathiterator(filenames)
        find_format = dxtbx.datablock.FormatChecker()

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
                if unhandled is not None:
                    unhandled.append(filename)
            elif format_class.is_abstract():
                logger.debug(
                    f"Image file {filename} appears to be a '{format_class.__name__}', but this is an abstract Format"
                )
                # Invalid format class found?
                if unhandled is not None:
                    unhandled.append(filename)
            elif issubclass(format_class, FormatMultiImage):
                imageset = format_class.get_imageset(
                    os.path.abspath(filename), format_kwargs=format_kwargs
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

        # Now, build experiments from these files. Duplicating the logic of
        # the previous implementation:
        # - FormatMultiImage files each have their own ImageSet
        # - Every set of images forming a scan goes into its own ImageSequence
        # - Any consecutive still frames that share any metadata with the
        #   previous still fram get collected into one ImageSet

        # Treat each format as a separate datablock
        for format_class, records in format_groups.items():
            if issubclass(format_class, FormatMultiImage):
                for imageset in records:
                    experiments.extend(
                        ExperimentListFactory.from_imageset_and_crystal(
                            imageset, crystal=None, load_models=load_models
                        )
                    )
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
            for imageset in imagesets:
                experiments.extend(
                    ExperimentListFactory.from_imageset_and_crystal(
                        imageset, crystal=None, load_models=load_models
                    )
                )

        # Assign indices to experiments
        for j, expt in enumerate(experiments):
            expt.index = j

        return experiments

    @staticmethod
    def from_imageset_and_crystal(imageset, crystal, load_models=True):
        """Load an experiment list from an imageset and crystal."""
        if isinstance(imageset, ImageSequence):
            return ExperimentListFactory.from_sequence_and_crystal(
                imageset, crystal, load_models
            )
        else:
            return ExperimentListFactory.from_stills_and_crystal(
                imageset, crystal, load_models
            )

    @staticmethod
    def from_sequence_and_crystal(imageset, crystal, load_models=True):
        """Create an experiment list from sequence and crystal."""

        assert isinstance(imageset, ImageSequence)

        experiments = ExperimentList()

        if load_models:
            # if imagesequence is still images, make one experiment for each
            # all referencing into the same image set
            if imageset.get_scan().is_still():
                start, end = imageset.get_scan().get_array_range()
                for j in range(start, end):
                    subset = imageset[j : j + 1]
                    experiments.append(
                        Experiment(
                            imageset=imageset,
                            beam=imageset.get_beam(),
                            detector=imageset.get_detector(),
                            goniometer=imageset.get_goniometer(),
                            scan=subset.get_scan(),
                            crystal=crystal,
                        )
                    )
            else:
                experiments.append(
                    Experiment(
                        imageset=imageset,
                        beam=imageset.get_beam(),
                        detector=imageset.get_detector(),
                        goniometer=imageset.get_goniometer(),
                        scan=imageset.get_scan(),
                        crystal=crystal,
                    )
                )

            return experiments

        else:
            return ExperimentList([Experiment(imageset=imageset, crystal=crystal)])

    @staticmethod
    def from_stills_and_crystal(imageset, crystal, load_models=True):
        """Create an experiment list from stills and crystal."""
        experiments = ExperimentList()
        if load_models:
            for i in range(len(imageset)):
                experiments.append(
                    Experiment(
                        imageset=imageset[i : i + 1],
                        beam=imageset.get_beam(i),
                        detector=imageset.get_detector(i),
                        goniometer=imageset.get_goniometer(i),
                        scan=imageset.get_scan(i),
                        crystal=crystal,
                    )
                )
        else:
            for i in range(len(imageset)):
                experiments.append(
                    Experiment(imageset=imageset[i : i + 1], crystal=crystal)
                )
        return experiments

    @staticmethod
    def from_datablock_and_crystal(datablock, crystal, load_models=True):
        """Load an experiment list from a datablock."""

        # Initialise the experiment list
        experiments = ExperimentList()

        # If we have a list, loop through
        if isinstance(datablock, list):
            for db in datablock:
                experiments.extend(
                    ExperimentListFactory.from_datablock_and_crystal(
                        db, crystal, load_models
                    )
                )
            return experiments

        # Add all the imagesets
        for imageset in datablock.extract_imagesets():
            experiments.extend(
                ExperimentListFactory.from_imageset_and_crystal(
                    imageset, crystal, load_models
                )
            )

        # Check the list is consistent
        assert experiments.is_consistent()

        return experiments

    @staticmethod
    def from_dict(obj, check_format=True, directory=None):
        """
        Load an experiment list from a dictionary.

        Args:
            obj (dict):
                Dictionary containing either ExperimentList or DataBlock
                structure.
            check_format (bool):
                If True, the file will be read to verify metadata.
            directory (str):

        Returns:
            ExperimentList: The dictionary converted
        """

        try:
            experiments = ExperimentList()
            for db in DataBlockFactory.from_dict(
                obj, check_format=check_format, directory=directory
            ):
                experiments.extend(
                    ExperimentListFactory.from_datablock_and_crystal(db, None)
                )
        except Exception:
            experiments = None

        # Decode the experiments from the dictionary
        if experiments is None:
            experiments = ExperimentListDict(
                obj, check_format=check_format, directory=directory
            ).decode()

        # Check the list is consistent
        assert experiments.is_consistent()

        # Assign indices to experiments
        for j, expt in enumerate(experiments):
            expt.index = j

        return experiments

    @staticmethod
    def from_json(text, check_format=True, directory=None):
        """Load an experiment list from JSON."""
        return ExperimentListFactory.from_dict(
            json.loads(text),
            check_format=check_format,
            directory=directory,
        )

    @staticmethod
    def from_json_file(filename, check_format=True):
        """Load an experiment list from a json file."""
        filename = os.path.abspath(filename)
        directory = os.path.dirname(filename)
        with open(filename, "r") as infile:
            return ExperimentListFactory.from_json(
                infile.read(), check_format=check_format, directory=directory
            )

    @staticmethod
    def from_pickle_file(filename):
        """Decode an experiment list from a pickle file."""
        with open(filename, "rb") as infile:
            obj = pickle.load(infile)
        assert isinstance(obj, ExperimentList)
        return obj

    @staticmethod
    def from_xds(xds_inp, xds_other):
        """Generate an experiment list from XDS files."""
        # Get the sequence from the XDS files
        sequence = xds.to_imageset(xds_inp, xds_other)

        # Get the crystal from the XDS files
        crystal = xds.to_crystal(xds_other)

        # Create the experiment list
        experiments = ExperimentListFactory.from_imageset_and_crystal(sequence, crystal)

        # Set the crystal in the experiment list
        assert len(experiments) == 1

        # Assign the index
        experiments[0].index = 0

        return experiments

    @staticmethod
    def from_serialized_format(filename, check_format=True):
        """Try to load the experiment list from a serialized format."""

        if hasattr(filename, "__fspath__"):
            filename = filename.__fspath__()  # unwrap PEP-519-style objects

        # First try as a JSON file
        try:
            return ExperimentListFactory.from_json_file(filename, check_format)
        except (FileNotFoundError, PermissionError):
            raise
        except Exception:
            pass

        # Now try as a pickle file
        return ExperimentListFactory.from_pickle_file(filename)

    @staticmethod
    def from_templates(templates, **kwargs):
        """Import an experiment list from templates"""
        importer = DataBlockTemplateImporter(templates, **kwargs)
        experiments = ExperimentList()
        for db in importer.datablocks:
            experiments.extend(
                ExperimentListFactory.from_datablock_and_crystal(db, None)
            )

        # Assign indices to experiments
        for j, expt in enumerate(experiments):
            expt.index = j

        return experiments


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


def _iterate_with_previous(iterable):
    """Convenience iterator to give pairs of (previous, next) items"""
    previous = None
    for val in iterable:
        yield (previous, val)
        previous = val


def _groupby_template_is_none(records):
    # type: (Iterable[ImageMetadataRecord]) -> Generator[List[ImageMetadataRecord]]
    """Specialization of groupby that groups records by format=None"""
    for _, group in itertools.groupby(
        enumerate(records), key=lambda x: -1 if x[1].template is None else x[0]
    ):
        yield list(x[1] for x in group)


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
        scan_tolerance: Fraction of oscillation range to tolerate
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
        x.filename if get_url_scheme(x.filename) else os.path.abspath(x.filename)
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


class ExperimentListTemplateImporter(object):
    """[DEPRECATED] Import an experiment list from a template.

    To be removed after DIALS v3.5 release branch is made; 3.5 v3.5.0"""

    def __init__(self, templates, **kwargs):
        warnings.warn(
            "ExperimentListTemplateImporter is deprecated; Please use ExperimentList.from_templates.",
            DeprecationWarning,
            stacklevel=2,
        )
        kwargs.pop("verbose", None)
        self.experiments = ExperimentListFactory.from_templates(templates, **kwargs)
