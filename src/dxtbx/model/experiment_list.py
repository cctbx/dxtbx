from __future__ import annotations

import collections
import copy
import errno
import importlib.metadata
import itertools
import json
import logging
import operator
import os
import pickle
import sys
from collections.abc import Callable, Generator, Iterable
from typing import Any

import natsort
from tqdm import tqdm

import dxtbx
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.format.image import ImageBool, ImageDouble
from dxtbx.format.Registry import get_format_class_for_file
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
from dxtbx.sequence_filenames import (
    locate_files_matching_template_string,
    template_image_range,
    template_regex,
    template_string_number_index,
)
from dxtbx.serialize import xds
from dxtbx.serialize.filename import resolve_path
from dxtbx.util import get_url_scheme

__all__ = [
    "BeamComparison",
    "DetectorComparison",
    "ExperimentListFactory",
    "GoniometerComparison",
]


logger = logging.getLogger(__name__)

# REMOVE and inline when Python 3.10 is minimum
if sys.version_info < (3, 10):
    scaling_model_entry_points = importlib.metadata.entry_points().get(
        "dxtbx.scaling_model_ext", []
    )
else:
    scaling_model_entry_points = importlib.metadata.entry_points(
        group="dxtbx.scaling_model_ext"
    )


class InvalidExperimentListError(RuntimeError):
    """
    Indicates an error whilst validating the experiment list.

    This means that there is some structural problem that prevents the given data
    from representing a well-formed experiment list. This doesn't indicate e.g.
    some problem with the data or model consistency.
    """


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


class ExperimentListDict:
    """A helper class for serializing the experiment list to dictionary (needed
    to save the experiment list to JSON format."""

    def __init__(self, obj, check_format=True, directory=None):
        """Initialise. Copy the dictionary."""
        # Basic check: This is a dict-like object. This can happen if e.g. we
        # were passed a DataBlock list instead of an ExperimentList dictionary
        if isinstance(obj, list) or not hasattr(obj, "get"):
            raise InvalidExperimentListError(f"Expected dictionary, not {type(obj)}")

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
        loaded and the experiment is rewritten with an integer pointing
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

    def _load_pickle_path(
        self, imageset_data: dict, param: str
    ) -> tuple[str | None, Any]:
        """
        Read a filename from an imageset dict and load if available. This is
        used to load mask, gain, pedestal and offset maps. In some situations
        (such as tests) these files are not available, in which case the
        filename is kept but the data is None.

        Args:
            imageset_data: The dictionary holding imageset information
            param: The key name to lookup in the imageset dictionary

        Returns:
            A tuple of (filename, data) where data has been loaded from
            the pickle file, or is None if the file is inaccessible. If there
            is no key entry then ("", None) is returned.
        """
        if param not in imageset_data:
            return "", None

        filename = resolve_path(imageset_data[param], directory=self._directory)
        data = None
        if filename:
            try:
                with open(filename, "rb") as fh:
                    data = pickle.load(fh, encoding="bytes")
            except OSError:
                pass
        else:
            filename = ""

        return filename, data

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

        # If dx, dy maps are expected then they must be loaded even when
        # self._check_format == False, because they affect the operation of
        # programs (dials.index, dials.refine) that do not need the image data.
        if (dx_filename or dy_filename) and not all((dx, dx)):
            raise RuntimeError(
                f"dx ({dx_filename}) and dy ({dy_filename}) maps are expected"
            )

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
                    assert abs(o[1] - (s[1])) < 1e-7
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
        for entry_point in scaling_model_entry_points:
            if entry_point.name == obj["__id__"]:
                return entry_point.load().from_dict(obj)


def _experimentlist_from_file(filename, directory=None):
    """Load a model dictionary from a file."""
    filename = resolve_path(filename, directory=directory)
    try:
        with open(filename) as infile:
            return json.load(infile)
    except OSError:
        raise OSError("unable to read file, %s" % filename)


class ExperimentListFactory:
    """A class to help instantiate experiment lists."""

    @staticmethod
    def from_args(
        args: list[str], unhandled: list[str] | None = None, check_format: bool = True
    ) -> ExperimentList:
        """Try to load serialised experiments from any recognised format."""

        # Create a list for unhandled arguments
        if unhandled is None:
            unhandled = []

        experiments = ExperimentList()

        # Try to load from serialized formats
        for filename in args:
            try:
                experiments.extend(
                    ExperimentListFactory.from_serialized_format(
                        filename, check_format=check_format
                    )
                )
                logger.debug(f"Loaded experiments from {filename}")
            except Exception as e:
                logger.debug(f"Could not load experiments from {filename}: {e}")
                unhandled.append(filename)
                raise

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
    ) -> ExperimentList:
        """Create a list of data blocks from a list of directory or file names."""
        experiments = ExperimentList()

        # Cast filenames to a list from whatever iterator they are
        filenames = list(filenames)

        # Process each file given by this path list
        to_process = _openingpathiterator(filenames)
        find_format = FormatChecker()

        format_groups = collections.defaultdict(list)
        if format_kwargs is None:
            format_kwargs = {}

        if os.isatty and len(filenames) > 1 and "DIALS_NOBANNER" not in os.environ:
            filename_iter = tqdm(to_process, total=len(filenames), file=sys.stdout)
        else:
            filename_iter = to_process

        for filename in filename_iter:
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
                format_groups[format_class].append(imageset)
                logger.debug("Loaded file: %s", filename)
            else:
                format_object = format_class(filename, **format_kwargs)
                meta = ImageMetadataRecord.from_format(format_object)
                assert meta.filename == filename

                # Add this entry to our table of formats
                format_groups[format_class].append(meta)
                logger.debug("Loaded metadata of file: %s", filename)

        # Now, build experiments from these files. Duplicating the logic of
        # the previous implementation:
        # - FormatMultiImage files each have their own ImageSet
        # - Every set of images forming a scan goes into its own ImageSequence
        # - Any consecutive still frames that share any metadata with the
        #   previous still fram get collected into one ImageSet

        all_tof = False
        for format_class, records in format_groups.items():
            for i in records:
                try:  # records can be ImageMetadataRecord or ImageSequence
                    scan = i.get_scan()
                    if scan is not None and scan.has_property("time_of_flight"):
                        all_tof = True
                    elif all_tof:
                        raise RuntimeError(
                            "Cannot process mix of ToF and non ToF experiments"
                        )
                except AttributeError:
                    if all_tof:
                        raise RuntimeError(
                            "Cannot process mix of ToF and non ToF experiments"
                        )

        # Treat each format as a separate block of data
        for format_class, records in format_groups.items():
            if issubclass(format_class, FormatMultiImage):
                if all_tof:
                    _merge_sequence_model_metadata(
                        records,
                        compare_beam=compare_beam,
                        compare_detector=compare_detector,
                        compare_goniometer=compare_goniometer,
                    )

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
            imagesets = list(
                _convert_to_imagesets(records, format_class, format_kwargs)
            )

            assert imagesets, "Got no imagesets when constructing ExperimentList?"
            for imageset in imagesets:
                experiments.extend(
                    ExperimentListFactory.from_imageset_and_crystal(
                        imageset, crystal=None, load_models=load_models
                    )
                )

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
                for j in range(len(imageset)):
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

        experiments = ExperimentListDict(
            obj, check_format=check_format, directory=directory
        ).decode()

        # Check the list is consistent
        assert experiments.is_consistent()

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
        try:
            with open(filename) as infile:
                return ExperimentListFactory.from_json(
                    infile.read(), check_format=check_format, directory=directory
                )
        except UnicodeDecodeError:
            raise InvalidExperimentListError(
                f"Cannot interpret {filename} as an ExperimentList"
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

        return experiments

    @staticmethod
    def from_serialized_format(filename, check_format=True):
        """Try to load the experiment list from a serialized format."""

        if hasattr(filename, "__fspath__"):
            filename = filename.__fspath__()  # unwrap PEP-519-style objects

        return ExperimentListFactory.from_json_file(filename, check_format)

    @staticmethod
    def from_templates(templates, **kwargs):
        """Import an experiment list from templates"""
        assert "verbose" not in kwargs, "The verbose parameter has been removed"
        assert len(templates) > 0

        experiments = ExperimentList()
        find_format = FormatChecker()

        # For each template do an import
        for template in templates:
            template = os.path.normpath(template)
            filenames = sorted(locate_files_matching_template_string(template))
            if len(filenames):
                logger.debug(
                    "The following files matched the template string:\n%s",
                    "\n".join(f" {p}" for p in filenames),
                )

            # Check if we've matched any filenames
            if len(filenames) == 0:
                raise ValueError(f"Template '{template}' does not match any files")

            # Get the format from the first image
            format_class = find_format.find_format(filenames[0])

            # Verify this makes sense
            if format_class is None:
                raise ValueError(f"Image file {filenames[0]} format is unknown")
            elif format_class.is_abstract():
                raise ValueError(
                    f"Image file {filenames[0]} appears to be a '{type(format_class).__name__}', but this is an abstract Format"
                )
            else:
                image_range = kwargs.get("image_range")
                if image_range:
                    first, last = image_range
                else:
                    first, last = template_image_range(template)

                if not kwargs.get("allow_incomplete_sequences", False):
                    if "#" in template:
                        # Check all images in range are present - if allowed
                        i0, i1 = template_string_number_index(template)
                        prefix = template[:i0]
                        suffix = template[i1:]
                        all_numbers = {
                            int(f.replace(prefix, "").replace(suffix, ""))
                            for f in filenames
                        }
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
                    else:
                        print(
                            "Warning: Using only one template file: %s. \n "
                            "`allow_incomplete_sequence` has no effect" % template
                        )

                # Read the image
                fmt = format_class(filenames[0], **(kwargs.get("format_kwargs", {})))

                # Update the image range
                image_range = (first, last)
                scan = fmt.get_scan()
                scan.set_image_range(image_range)

                # Create the sequence and experiment
                imageset = dxtbx.imageset.ImageSetFactory.make_sequence(
                    template,
                    list(range(first, last + 1)),
                    format_class,
                    fmt.get_beam(),
                    fmt.get_detector(),
                    fmt.get_goniometer(),
                    scan,
                    format_kwargs=kwargs.get("format_kwargs"),
                )
                experiments.extend(
                    ExperimentListFactory.from_imageset_and_crystal(
                        imageset,
                        crystal=None,
                        load_models=True,
                    )
                )
        return experiments


class ImageMetadataRecord:
    """Object to store metadata information.

    This is used whilst building the experiment lists. The metadata for each
    image can be read once, and then any grouping/deduplication can happen
    later, without re-opening the original file.
    """

    def __init__(
        self,
        beam: dxtbx.model.Beam | None = None,
        detector: dxtbx.model.Detector | None = None,
        goniometer: dxtbx.model.Goniometer | None = None,
        scan: dxtbx.model.Scan | None = None,
        template: str | None = None,
        filename: str | None = None,
        index: int | None = None,
    ):
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
        other_record: ImageMetadataRecord,
        compare_beam: Callable = operator.__eq__,
        compare_detector: Callable = operator.__eq__,
        compare_goniometer: Callable = operator.__eq__,
    ) -> bool:
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
    def from_format(cls, fmt: Format) -> Any:
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


def _groupby_template_is_none(
    records: Iterable[ImageMetadataRecord],
) -> Generator[list[ImageMetadataRecord], None, None]:
    """Specialization of groupby that groups records by format=None"""
    for _, group in itertools.groupby(
        enumerate(records), key=lambda x: -1 if x[1].template is None else x[0]
    ):
        yield [x[1] for x in group]


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
    paths = collections.deque((True, x) for x in natsort.natsorted(pathnames))

    while paths:
        # Get the next path from the queue
        (do_recurse, pathname) = paths.popleft()
        pathname = os.fspath(pathname)
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
    records: Iterable[ImageMetadataRecord],
    compare_beam: Callable | None = None,
    compare_detector: Callable | None = None,
    compare_goniometer: Callable | None = None,
):
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


def _merge_sequence_model_metadata(
    records: Iterable[ImageSequence],
    compare_beam: Callable | None = None,
    compare_detector: Callable | None = None,
    compare_goniometer: Callable | None = None,
):
    record_altered = False
    for prev, record in _iterate_with_previous(records):
        if prev is None:
            continue

        record_altered = False
        record_beam = record.get_beam()
        record_detector = record.get_detector()
        record_goniometer = record.get_goniometer()
        prev_beam = prev.get_beam()
        prev_detector = prev.get_detector()
        prev_goniometer = prev.get_goniometer()

        if record_beam is not prev_beam and compare_beam(record_beam, prev_beam):
            record.set_beam(prev_beam)
            record_altered = True
        if record_detector is not prev_detector and compare_detector(
            record_detector, prev_detector
        ):
            record.set_detector(prev_detector)
            record_altered = True
        if record_goniometer is not prev_goniometer and compare_goniometer(
            record_goniometer, prev_goniometer
        ):
            record.set_goniometer(prev_goniometer)
            record_altered = True

    return record_altered


def _merge_scans(
    records: Iterable[ImageMetadataRecord], scan_tolerance: float | None = None
) -> list[ImageMetadataRecord]:
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
                logger.debug(
                    "  Failed to merge record %s with previous - writing new scan",
                    str(e),
                )
            else:
                # If we appended, then we don't need to keep this record's scan
                record.scan = prev.scan
                logger.debug("  Appended record %s to previous", record)
                continue
        merged_records.append(record)

    logger.debug("Result of merging record scans: %d records", len(merged_records))
    return merged_records


def _convert_to_imagesets(
    records: Iterable[ImageMetadataRecord],
    format_class: type[Format],
    format_kwargs: dict | None = None,
) -> Generator[dxtbx.imageset.ImageSet, None, None]:
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


def _create_imageset(
    records: Iterable[ImageMetadataRecord],
    format_class: type[Format],
    format_kwargs: dict | None = None,
) -> dxtbx.imageset.ImageSet:
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
    # Everything should have a filename
    assert all(x.filename for x in records)
    # Extract the filenames from the records
    filenames = [
        x.filename if get_url_scheme(x.filename) else os.path.abspath(x.filename)
        for x in records
        if x.filename
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


def _create_imagesequence(
    record: ImageMetadataRecord,
    format_class: type[Format],
    format_kwargs: dict | None = None,
) -> dxtbx.imageset.ImageSequence:
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
    assert record.scan
    assert record.template
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
