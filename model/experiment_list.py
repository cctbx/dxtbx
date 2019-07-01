from __future__ import absolute_import, division, print_function

from builtins import range
import copy
import json
import os
import pkg_resources
import warnings

from dxtbx.datablock import (
    BeamComparison,
    DataBlockFactory,
    DataBlockTemplateImporter,
    DetectorComparison,
    GoniometerComparison,
    SweepDiff,
)
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.format.image import ImageBool, ImageDouble
from dxtbx.imageset import ImageGrid, ImageSet, ImageSetFactory, ImageSweep
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
from dxtbx.serialize import xds
from dxtbx.serialize.filename import resolve_path
from dxtbx.serialize.load import _decode_dict
from dxtbx.sweep_filenames import template_image_range
import six.moves.cPickle as pickle

try:
    from typing import Any, Dict, Optional, Tuple
except ImportError:
    pass

__all__ = [
    "BeamComparison",
    "DetectorComparison",
    "ExperimentListFactory",
    "GoniometerComparison",
    "SweepDiff",
]


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
        """ Initialise. Copy the dictionary. """
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
        """ Helper function. Extract the models. """
        """if name == imageset: Extract imageset objects from the source.

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

        # Return the model list
        return mlist

    def _load_pickle_path(self, imageset_data, param):
        # type: (Dict, str) -> Tuple[Optional[str], Any]
        """Read a filename from an imageset dict and load if required.

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
                return filename, pickle.load(fh)
        return filename or "", None

    def decode(self):
        """ Decode the dictionary into a list of experiments. """
        # Extract all the experiments

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

            key = (eobj.get("imageset"), eobj.get("scan"))

            imageset = None
            try:
                imageset = imagesets[key]  # type: ImageSet
            except KeyError:
                # This imageset hasn't been loaded yet - create it
                imageset_data = self._lookup_model("imageset", eobj)

                # Create the imageset from the input data
                if imageset_data is not None:
                    if "params" in imageset_data:
                        format_kwargs = imageset_data["params"]
                    else:
                        format_kwargs = {}

                    # Load the external lookup data
                    mask_filename, mask = self._load_pickle_path(imageset_data, "mask")
                    gain_filename, gain = self._load_pickle_path(imageset_data, "gain")
                    pedestal_filename, pedestal = self._load_pickle_path(
                        imageset_data, "pedestal"
                    )
                    dx_filename, dx = self._load_pickle_path(imageset_data, "dx")
                    dy_filename, dy = self._load_pickle_path(imageset_data, "dy")

                    if imageset_data["__id__"] == "ImageSet":
                        imageset = self._make_stills(
                            imageset_data, format_kwargs=format_kwargs
                        )
                    elif imageset_data["__id__"] == "ImageGrid":
                        imageset = self._make_grid(
                            imageset_data, format_kwargs=format_kwargs
                        )
                    elif imageset_data["__id__"] == "ImageSweep":
                        imageset = self._make_sweep(
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
                        if isinstance(imageset, ImageSweep):
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

                # Add the imageset to the dict - even if empty - as this will
                # prevent a duplicated attempt at reconstruction
                imagesets[key] = imageset

            # Append the experiment
            el.append(
                Experiment(
                    imageset=imageset,
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

        # Return the experiment list
        return el

    def _make_mem_imageset(self, imageset):
        """ Can't make a mem imageset from dict. """
        return None

    def _make_stills(self, imageset, format_kwargs=None):
        """ Make a still imageset. """
        filenames = [
            resolve_path(p, directory=self._directory) for p in imageset["images"]
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
        """ Make a still imageset. """
        grid_size = imageset["grid_size"]
        return ImageGrid.from_imageset(
            self._make_stills(imageset, format_kwargs=format_kwargs), grid_size
        )

    def _make_sweep(
        self,
        imageset,
        beam=None,
        detector=None,
        goniometer=None,
        scan=None,
        format_kwargs=None,
    ):
        """ Make an image sweep. """
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

        # Make a sweep from the input data
        return ImageSetFactory.make_sweep(
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
        """ Get the scaling model from a dictionary. """
        for entry_point in pkg_resources.iter_entry_points("dxtbx.scaling_model_ext"):
            if entry_point.name == obj["__id__"]:
                return entry_point.load().from_dict(obj)


def _experimentlist_from_file(filename, directory=None):
    """ Load a model dictionary from a file. """
    filename = resolve_path(filename, directory=directory)
    try:
        with open(filename, "r") as infile:
            return json.load(infile, object_hook=_decode_dict)
    except IOError:
        raise IOError("unable to read file, %s" % filename)


class ExperimentListDumper(object):
    """ A class to help writing JSON files. """

    def __init__(self, experiment_list):
        """ Initialise """
        warnings.warn(
            "class ExperimentListDumper() is deprecated. "
            "Use experiment_list.as_json(), experiment_list.as_pickle(), experiment_list.as_file() directly",
            DeprecationWarning,
            stacklevel=2,
        )
        assert experiment_list
        self.as_json = experiment_list.as_json
        self.as_pickle = experiment_list.as_pickle
        self.as_file = experiment_list.as_file


class ExperimentListFactory(object):
    """ A class to help instantiate experiment lists. """

    @staticmethod
    def from_args(args, verbose=False, unhandled=None):
        """ Try to load experiment from any recognised format. """

        # Create a list for unhandled arguments
        if unhandled is None:
            unhandled = []

        experiments = ExperimentList()
        ## First try as image files
        # experiments = ExperimentListFactory.from_datablock(
        # DataBlockFactory.from_args(args, verbose, unhandled1))

        # Try to load from serialized formats
        for filename in args:
            try:
                experiments.extend(
                    ExperimentListFactory.from_serialized_format(filename)
                )
                if verbose:
                    print("Loaded experiments from %s" % filename)
            except Exception as e:
                if verbose:
                    print("Could not load experiments from %s: %s" % (filename, str(e)))
                unhandled.append(filename)

        # Return the experiments
        return experiments

    @staticmethod
    def from_filenames(
        filenames,
        verbose=False,
        unhandled=None,
        compare_beam=None,
        compare_detector=None,
        compare_goniometer=None,
        scan_tolerance=None,
        format_kwargs=None,
        load_models=True,
    ):
        """ Create a list of data blocks from a list of directory or file names. """
        experiments = ExperimentList()
        for db in DataBlockFactory.from_filenames(
            filenames,
            verbose=verbose,
            unhandled=unhandled,
            compare_beam=compare_beam,
            compare_detector=compare_detector,
            compare_goniometer=compare_goniometer,
            scan_tolerance=scan_tolerance,
            format_kwargs=format_kwargs,
        ):
            experiments.extend(
                ExperimentListFactory.from_datablock_and_crystal(db, None, load_models)
            )
        return experiments

    @staticmethod
    def from_imageset_and_crystal(imageset, crystal, load_models=True):
        """ Load an experiment list from an imageset and crystal. """
        if isinstance(imageset, ImageSweep):
            return ExperimentListFactory.from_sweep_and_crystal(
                imageset, crystal, load_models
            )
        else:
            return ExperimentListFactory.from_stills_and_crystal(
                imageset, crystal, load_models
            )

    @staticmethod
    def from_sweep_and_crystal(imageset, crystal, load_models=True):
        """ Create an experiment list from sweep and crystal. """
        if load_models:
            return ExperimentList(
                [
                    Experiment(
                        imageset=imageset,
                        beam=imageset.get_beam(),
                        detector=imageset.get_detector(),
                        goniometer=imageset.get_goniometer(),
                        scan=imageset.get_scan(),
                        crystal=crystal,
                    )
                ]
            )
        else:
            return ExperimentList([Experiment(imageset=imageset, crystal=crystal)])

    @staticmethod
    def from_stills_and_crystal(imageset, crystal, load_models=True):
        """ Create an experiment list from stills and crystal. """
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
        """ Load an experiment list from a datablock. """

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

        # Return the experiments
        return experiments

    @staticmethod
    def from_dict(obj, check_format=True, directory=None):
        """Load an experiment list from a dictionary.

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

        # Return the experiments
        return experiments

    @staticmethod
    def from_json(text, check_format=True, directory=None):
        """ Load an experiment list from JSON. """
        return ExperimentListFactory.from_dict(
            json.loads(text, object_hook=_decode_dict),
            check_format=check_format,
            directory=directory,
        )

    @staticmethod
    def from_json_file(filename, check_format=True):
        """ Load an experiment list from a json file. """
        filename = os.path.abspath(filename)
        directory = os.path.dirname(filename)
        with open(filename, "r") as infile:
            return ExperimentListFactory.from_json(
                infile.read(), check_format=check_format, directory=directory
            )

    @staticmethod
    def from_pickle_file(filename):
        """ Decode an experiment list from a pickle file. """
        with open(filename, "rb") as infile:
            obj = pickle.load(infile)
        assert isinstance(obj, ExperimentList)
        return obj

    @staticmethod
    def from_xds(xds_inp, xds_other):
        """ Generate an experiment list from XDS files. """
        # Get the sweep from the XDS files
        sweep = xds.to_imageset(xds_inp, xds_other)

        # Get the crystal from the XDS files
        crystal = xds.to_crystal(xds_other)

        # Create the experiment list
        experiments = ExperimentListFactory.from_imageset_and_crystal(sweep, crystal)

        # Set the crystal in the experiment list
        assert len(experiments) == 1

        # Return the experiment list
        return experiments

    @staticmethod
    def from_serialized_format(filename, check_format=True):
        """ Try to load the experiment list from a serialized format. """

        # First try as a JSON file
        try:
            return ExperimentListFactory.from_json_file(filename, check_format)
        except Exception:
            pass

        # Now try as a pickle file
        return ExperimentListFactory.from_pickle_file(filename)


class ExperimentListTemplateImporter(object):
    """ A class to import an experiment list from a template. """

    def __init__(self, templates, verbose=False, **kwargs):
        importer = DataBlockTemplateImporter(templates, verbose=verbose, **kwargs)
        self.experiments = ExperimentList()
        for db in importer.datablocks:
            self.experiments.extend(
                ExperimentListFactory.from_datablock_and_crystal(db, None)
            )
