from __future__ import absolute_import, division, print_function

import collections
import json
from copy import deepcopy
from os.path import abspath, dirname, splitext

import pkg_resources
import six.moves.cPickle as pickle

from dxtbx.datablock import (
    AutoEncoder,
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

        self._obj = deepcopy(obj)
        self._check_format = check_format
        self._directory = directory

    def decode(self):
        """ Decode the dictionary into a list of experiments. """

        # If this doesn't claim to be an ExperimentList, don't even try
        if not self._obj.get("__id__", None) == "ExperimentList":
            raise InvalidExperimentListError(
                "Expected __id__ 'ExperimentList', but found {}".format(
                    repr(self._obj.get("__id__"))
                )
            )

        # Extract lists of models referenced by experiments
        self._blist = self._extract_models("beam")
        self._dlist = self._extract_models("detector")
        self._glist = self._extract_models("goniometer")
        self._slist = self._extract_models("scan")
        self._clist = self._extract_models("crystal")
        self._plist = self._extract_models("profile")
        self._scalelist = self._extract_models("scaling_model")

        # Go through all the imagesets and make sure the dictionary
        # references by an index rather than a file path.
        self._ilist = self._extract_imagesets()

        # Extract all the experiments
        return self._extract_experiments()

    def _extract_models(self, name):
        """ Helper function. Extract the models. """

        # The from dict function
        from_dict = getattr(self, "_%s_from_dict" % name)

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
                        from_dict(ExperimentListDict._from_file(value, self._directory))
                    )
                eobj[name] = mmap[value]
            elif not isinstance(value, int):
                raise TypeError("expected int or str, got %s" % type(value))

        # Return the model list
        return mlist

    def _extract_imagesets(self):
        """Extract imageset objects from the source.

        This function does resolving of an (old) method of imageset lookup
        e.g. it was valid to have a string as the imageset value in an
        experiment instead of an int - in which case the imageset was
        loaded from the named file in the target directory.

        If any experiments point to a file in this way, the imageset is
        loaded and the experiment is rewritted with an integer pointing
        to the new ImageSet in the returned list.

        Returns:
            List[ImageSet]:
                The ordered list of ImageSets that the Experiment index
                points to.
        """

        # Extract all the model list
        mlist = self._obj.get("imageset", [])

        # Dictionaries for file mappings
        mmap = {}

        # For each experiment, check the imageset is not specified by
        # a path, if it is then get the dictionary of the imageset
        # and insert it into the list. Replace the path reference
        # with an index
        for eobj in self._obj["experiment"]:
            value = eobj.get("imageset")
            if value is None:
                continue
            elif isinstance(value, str):
                if value not in mmap:
                    mmap[value] = len(mlist)
                    mlist.append(ExperimentListDict._from_file(value, self._directory))
                eobj["imageset"] = mmap[value]
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
            return None, None

        filename = resolve_path(imageset_data[param], directory=self._directory)
        if self._check_format and filename:
            with open(filename, "rb") as fh:
                return filename, pickle.load(fh)
        return filename, None

    def _extract_experiments(self):
        """ Helper function. Extract the experiments. """
        # Map of imageset/scan pairs
        imagesets = {}

        # For every experiment, use the given input to create
        # a sensible experiment.
        el = ExperimentList()
        for eobj in self._obj["experiment"]:

            # Get the models
            identifier = eobj.get("identifier", "")
            beam = ExperimentListDict.model_or_none(self._blist, eobj, "beam")
            detector = ExperimentListDict.model_or_none(self._dlist, eobj, "detector")
            goniometer = ExperimentListDict.model_or_none(
                self._glist, eobj, "goniometer"
            )
            scan = ExperimentListDict.model_or_none(self._slist, eobj, "scan")
            crystal = ExperimentListDict.model_or_none(self._clist, eobj, "crystal")
            profile = ExperimentListDict.model_or_none(self._plist, eobj, "profile")
            scaling_model = ExperimentListDict.model_or_none(
                self._scalelist, eobj, "scaling_model"
            )
            key = (eobj.get("imageset"), eobj.get("scan"))

            imageset = None
            try:
                imageset = imagesets[key]  # type: ImageSet
            except KeyError:
                # This imageset hasn't been loaded yet - create it
                imageset_data = ExperimentListDict.model_or_none(
                    self._ilist, eobj, "imageset"
                )  # type: dict

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
                        if mask_filename is None:
                            mask_filename = ""
                        if gain_filename is None:
                            gain_filename = ""
                        if pedestal_filename is None:
                            pedestal_filename = ""
                        if dx_filename is None:
                            dx_filename = ""
                        if dy_filename is None:
                            dy_filename = ""
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
                        elif isinstance(imageset, ImageSet):
                            for i in range(len(imageset)):
                                imageset.set_beam(beam, i)
                                imageset.set_detector(detector, i)
                                imageset.set_goniometer(goniometer, i)
                                imageset.set_scan(scan, i)
                        elif isinstance(imageset, ImageGrid):
                            for i in range(len(imageset)):
                                imageset.set_beam(beam, i)
                                imageset.set_detector(detector, i)
                                imageset.set_goniometer(goniometer, i)
                                imageset.set_scan(scan, i)
                        else:
                            pass

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

    @staticmethod
    def model_or_none(mlist, eobj, name):
        """ Get a model or None. """
        index = eobj.get(name)
        if index is not None:
            return mlist[index]
        return None

    @staticmethod
    def _beam_from_dict(obj):
        """ Get a beam from a dictionary. """
        return BeamFactory.from_dict(obj)

    @staticmethod
    def _detector_from_dict(obj):
        """ Get the detector from a dictionary. """
        return DetectorFactory.from_dict(obj)

    @staticmethod
    def _goniometer_from_dict(obj):
        """ Get the goniometer from a dictionary. """
        return GoniometerFactory.from_dict(obj)

    @staticmethod
    def _scan_from_dict(obj):
        """ Get the scan from a dictionary. """
        return ScanFactory.from_dict(obj)

    @staticmethod
    def _crystal_from_dict(obj):
        """ Get the crystal from a dictionary. """
        return CrystalFactory.from_dict(obj)

    @staticmethod
    def _profile_from_dict(obj):
        """ Get the profile from a dictionary. """
        return ProfileModelFactory.from_dict(obj)

    @staticmethod
    def _scaling_model_from_dict(obj):
        """ Get the scaling model from a dictionary. """
        for entry_point in pkg_resources.iter_entry_points("dxtbx.scaling_model_ext"):
            if entry_point.name == obj["__id__"]:
                return entry_point.load().from_dict(obj)

    @staticmethod
    def _from_file(filename, directory=None):
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
        self._experiment_list = experiment_list

    def as_json(self, filename=None, compact=False, split=False):
        """ Dump experiment list as json """
        # Get the dictionary and get the JSON string
        dictionary = self._experiment_list.to_dict()

        # Split into separate files
        if filename is not None and split:

            # Get lists of models by filename
            basepath = splitext(filename)[0]
            ilist = [
                ("%s_imageset_%d.json" % (basepath, i), d)
                for i, d in enumerate(dictionary["imageset"])
            ]
            blist = [
                ("%s_beam_%d.json" % (basepath, i), d)
                for i, d in enumerate(dictionary["beam"])
            ]
            dlist = [
                ("%s_detector_%d.json" % (basepath, i), d)
                for i, d in enumerate(dictionary["detector"])
            ]
            glist = [
                ("%s_goniometer_%d.json" % (basepath, i), d)
                for i, d in enumerate(dictionary["goniometer"])
            ]
            slist = [
                ("%s_scan_%d.json" % (basepath, i), d)
                for i, d in enumerate(dictionary["scan"])
            ]
            clist = [
                ("%s_crystal_%d.json" % (basepath, i), d)
                for i, d in enumerate(dictionary["crystal"])
            ]
            plist = [
                ("%s_profile_%d.json" % (basepath, i), d)
                for i, d in enumerate(dictionary["profile"])
            ]
            scalelist = [
                ("%s_scaling_model_%d.json" % (basepath, i), d)
                for i, d in enumerate(dictionary["scaling_model"])
            ]

            # Get the list of experiments
            edict = collections.OrderedDict(
                [("__id__", "ExperimentList"), ("experiment", dictionary["experiment"])]
            )

            # Set paths rather than indices
            for e in edict["experiment"]:
                if "imageset" in e:
                    e["imageset"] = ilist[e["imageset"]][0]
                if "beam" in e:
                    e["beam"] = blist[e["beam"]][0]
                if "detector" in e:
                    e["detector"] = dlist[e["detector"]][0]
                if "goniometer" in e:
                    e["goniometer"] = glist[e["goniometer"]][0]
                if "scan" in e:
                    e["scan"] = slist[e["scan"]][0]
                if "crystal" in e:
                    e["crystal"] = clist[e["crystal"]][0]
                if "profile" in e:
                    e["profile"] = plist[e["profile"]][0]
                if "scaling_model" in e:
                    e["scaling_model"] = scalelist[e["scaling_model"]][0]

            to_write = (
                ilist
                + blist
                + dlist
                + glist
                + slist
                + clist
                + plist
                + scalelist
                + [(filename, edict)]
            )
        else:
            to_write = [(filename, dictionary)]

        for fname, obj in to_write:
            if compact:
                separators = (",", ":")
                indent = None
            else:
                separators = None
                indent = 2
            text = json.dumps(
                obj,
                separators=separators,
                indent=indent,
                ensure_ascii=True,
                cls=AutoEncoder,
            )

            # If a filename is set then dump to file otherwise return string
            if fname is not None:
                with open(fname, "w") as outfile:
                    outfile.write(text)
            else:
                return text

    def as_pickle(self, filename=None, **kwargs):
        """ Dump experiment list as pickle. """
        # Get the pickle string
        text = pickle.dumps(self._experiment_list, protocol=pickle.HIGHEST_PROTOCOL)

        # Write the file
        if filename:
            with open(filename, "wb") as outfile:
                outfile.write(text)
        else:
            return text

    def as_file(self, filename, **kwargs):
        """ Dump experiment list as file. """
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
        filename = abspath(filename)
        directory = dirname(filename)
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
