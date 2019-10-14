from __future__ import absolute_import, division, print_function

import collections
import json
import os
import sys
import warnings
from builtins import range

import boost.python
import cctbx.crystal
from libtbx.containers import OrderedSet
from libtbx.utils import format_float_with_standard_uncertainty
from scitbx import matrix
from scitbx.array_family import flex

import six.moves.cPickle as pickle
from dxtbx.imageset import ImageGrid, ImageSequence, ImageSet
from dxtbx.model.beam import BeamFactory
from dxtbx.model.crystal import CrystalFactory
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model.profile import ProfileModelFactory
from dxtbx.model.scan import ScanFactory
from dxtbx_model_ext import (
    Beam,
    BeamBase,
    Crystal,
    CrystalBase,
    Detector,
    DetectorNode,
    Experiment,
    ExperimentList,
    Goniometer,
    GoniometerBase,
    KappaDirection,
    KappaGoniometer,
    KappaScanAxis,
    MosaicCrystalKabsch2010,
    MosaicCrystalSauter2014,
    MultiAxisGoniometer,
    OffsetParallaxCorrectedPxMmStrategy,
    OffsetPxMmStrategy,
    Panel,
    ParallaxCorrectedPxMmStrategy,
    PxMmStrategy,
    Scan,
    ScanBase,
    SimplePxMmStrategy,
    VirtualPanel,
    VirtualPanelFrame,
    get_mod2pi_angles_in_range,
    get_range_of_mod2pi_angles,
    is_angle_in_range,
    parallax_correction,
    parallax_correction_inv,
)

__all__ = (
    "Beam",
    "BeamBase",
    "BeamFactory",
    "Crystal",
    "CrystalBase",
    "CrystalFactory",
    "Detector",
    "DetectorFactory",
    "DetectorNode",
    "Experiment",
    "ExperimentList",
    "Goniometer",
    "GoniometerBase",
    "GoniometerFactory",
    "ImageGrid",
    "ImageSet",
    "ImageSequence",
    "KappaDirection",
    "KappaGoniometer",
    "KappaScanAxis",
    "MosaicCrystalKabsch2010",
    "MosaicCrystalSauter2014",
    "MultiAxisGoniometer",
    "OffsetParallaxCorrectedPxMmStrategy",
    "OffsetPxMmStrategy",
    "Panel",
    "ParallaxCorrectedPxMmStrategy",
    "ProfileModelFactory",
    "PxMmStrategy",
    "Scan",
    "ScanBase",
    "ScanFactory",
    "SimplePxMmStrategy",
    "VirtualPanel",
    "VirtualPanelFrame",
    "get_mod2pi_angles_in_range",
    "get_range_of_mod2pi_angles",
    "is_angle_in_range",
    "parallax_correction",
    "parallax_correction_inv",
)


@boost.python.inject_into(Detector)
class _(object):
    def iter_panels(self):
        """ Iterate through just the panels depth-first. """
        for obj in self.iter_preorder():
            if obj.is_panel():
                yield obj

    def iter_preorder(self):
        """ Iterate through the groups and panels depth-first. """
        stack = [self.hierarchy()]
        while stack:
            node = stack.pop()
            yield node
            if node.is_group():
                stack.extend(reversed(node))

    def iter_levelorder(self):
        """ Iterate through the groups and panels breadth-first. """
        queue = [self.hierarchy()]
        while queue:
            node = queue.pop(0)
            yield node
            if node.is_group():
                queue.extend(node)


@boost.python.inject_into(Crystal)
class _(object):
    def show(self, show_scan_varying=False, out=None):
        if out is None:
            out = sys.stdout
        print(self.as_str(show_scan_varying=show_scan_varying), file=out)

    def get_crystal_symmetry(self, assert_is_compatible_unit_cell=True):
        return cctbx.crystal.symmetry(
            unit_cell=self.get_unit_cell(),
            space_group=self.get_space_group(),
            assert_is_compatible_unit_cell=assert_is_compatible_unit_cell,
        )

    def as_str(self, show_scan_varying=False):
        uc = self.get_unit_cell().parameters()
        uc_sd = self.get_cell_parameter_sd()
        sg = str(self.get_space_group().info())
        umat = (
            matrix.sqr(self.get_U())
            .mathematica_form(format="% 5.4f", one_row_per_line=True)
            .splitlines()
        )
        bmat = (
            matrix.sqr(self.get_B())
            .mathematica_form(format="% 5.4f", one_row_per_line=True)
            .splitlines()
        )
        amat = (
            (matrix.sqr(self.get_U()) * matrix.sqr(self.get_B()))
            .mathematica_form(format="% 5.4f", one_row_per_line=True)
            .splitlines()
        )

        msg = ["Crystal:"]

        if len(uc_sd) != 0:
            cell_str = [
                format_float_with_standard_uncertainty(v, e, minimum=1.0e-5)
                for (v, e) in zip(uc, uc_sd)
            ]
            msg.append("    Unit cell: (" + ", ".join(cell_str) + ")")
        else:
            msg.append(
                "    Unit cell: " + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc
            )
        msg.append("    Space group: " + sg)
        msg.append("    U matrix:  " + umat[0])
        msg.append("               " + umat[1])
        msg.append("               " + umat[2])
        msg.append("    B matrix:  " + bmat[0])
        msg.append("               " + bmat[1])
        msg.append("               " + bmat[2])
        msg.append("    A = UB:    " + amat[0])
        msg.append("               " + amat[1])
        msg.append("               " + amat[2])
        if self.num_scan_points > 0:
            msg.append("    A sampled at " + str(self.num_scan_points) + " scan points")
            if show_scan_varying:
                for i in range(self.num_scan_points):
                    A = matrix.sqr(self.get_A_at_scan_point(i))
                    B = matrix.sqr(self.get_B_at_scan_point(i))
                    U = matrix.sqr(self.get_U_at_scan_point(i))
                    uc = self.get_unit_cell_at_scan_point(i).parameters()
                    umat = U.mathematica_form(
                        format="% 5.4f", one_row_per_line=True
                    ).splitlines()
                    bmat = B.mathematica_form(
                        format="% 5.4f", one_row_per_line=True
                    ).splitlines()
                    amat = A.mathematica_form(
                        format="% 5.4f", one_row_per_line=True
                    ).splitlines()
                    msg.append("  Scan point #%i:" % (i + 1))
                    msg.append(
                        "    Unit cell: "
                        + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc
                    )
                    msg.append("    U matrix:  " + umat[0])
                    msg.append("               " + umat[1])
                    msg.append("               " + umat[2])
                    msg.append("    B matrix:  " + bmat[0])
                    msg.append("               " + bmat[1])
                    msg.append("               " + bmat[2])
                    msg.append("    A = UB:    " + amat[0])
                    msg.append("               " + amat[1])
                    msg.append("               " + amat[2])
        return "\n".join(msg)

    def __str__(self):
        return self.as_str()

    def to_dict(self):
        """Convert the crystal model to a dictionary

        Returns:
            A dictionary of the parameters
        """
        # Get the real space vectors
        A = matrix.sqr(self.get_A()).inverse()
        real_space_a = (A[0], A[1], A[2])
        real_space_b = (A[3], A[4], A[5])
        real_space_c = (A[6], A[7], A[8])

        # Get the space group Hall symbol
        hall = self.get_space_group().info().type().hall_symbol()

        # Isoforms used for stills
        try:
            identified_isoform = self.identified_isoform
        except AttributeError:
            identified_isoform = None

        # Collect the information as a python dictionary
        xl_dict = collections.OrderedDict(
            [
                ("__id__", "crystal"),
                ("real_space_a", real_space_a),
                ("real_space_b", real_space_b),
                ("real_space_c", real_space_c),
                ("space_group_hall_symbol", hall),
            ]
        )

        if identified_isoform is not None:
            xl_dict["identified_isoform"] = identified_isoform

        # Add in scan points if present
        if self.num_scan_points > 0:
            A_at_scan_points = tuple(
                [self.get_A_at_scan_point(i) for i in range(self.num_scan_points)]
            )
            xl_dict["A_at_scan_points"] = A_at_scan_points

        # Add in covariance of B if present
        cov_B = tuple(self.get_B_covariance())
        if len(cov_B) != 0:
            xl_dict["B_covariance"] = cov_B

        # Add in covariance of B at scan points if present
        if self.num_scan_points > 0:
            try:
                cov_B_at_scan_points = tuple(
                    [
                        tuple(self.get_B_covariance_at_scan_point(i))
                        for i in range(self.num_scan_points)
                    ]
                )
                xl_dict["B_covariance_at_scan_points"] = cov_B_at_scan_points
            except RuntimeError:
                pass

        return xl_dict

    @staticmethod
    def from_dict(d):
        """Convert the dictionary to a crystal model

        Params:
            d The dictionary of parameters

        Returns:
            The crystal model
        """
        # If None, return None
        if d is None:
            return None

        # Check the version and id
        if str(d["__id__"]) != "crystal":
            raise ValueError('"__id__" does not equal "crystal"')

        # Extract from the dictionary
        real_space_a = d["real_space_a"]
        real_space_b = d["real_space_b"]
        real_space_c = d["real_space_c"]
        # str required to force unicode to ascii conversion
        space_group = str("Hall:" + d["space_group_hall_symbol"])
        xl = Crystal(
            real_space_a, real_space_b, real_space_c, space_group_symbol=space_group
        )

        # Isoforms used for stills
        try:
            xl.identified_isoform = d["identified_isoform"]
        except KeyError:
            pass

        # Extract scan point setting matrices, if present
        try:
            A_at_scan_points = d["A_at_scan_points"]
            xl.set_A_at_scan_points(A_at_scan_points)
        except KeyError:
            pass

        # Extract covariance of B, if present
        try:
            cov_B = d["B_covariance"]
            xl.set_B_covariance(cov_B)
        except KeyError:
            pass

        # Extract covariance of B at scan points, if present
        cov_B_at_scan_points = d.get("B_covariance_at_scan_points")
        if cov_B_at_scan_points is not None:
            cov_B_at_scan_points = flex.double(cov_B_at_scan_points).as_1d()
            cov_B_at_scan_points.reshape(flex.grid(xl.num_scan_points, 9, 9))
            xl.set_B_covariance_at_scan_points(cov_B_at_scan_points)

        return xl


@boost.python.inject_into(MosaicCrystalKabsch2010)
class _(object):
    def as_str(self, show_scan_varying=False):
        return "\n".join(
            (
                super(MosaicCrystalKabsch2010, self).as_str(
                    show_scan_varying=show_scan_varying
                ),
                "    Mosaicity:  %.6f" % self.get_mosaicity(),
            )
        )

    def to_dict(self):
        """Convert the crystal model to a dictionary

        Returns:
            A dictionary of the parameters

        """
        xl_dict = super(MosaicCrystalKabsch2010, self).to_dict()

        # Get the mosaicity
        mosaicity = self.get_mosaicity()
        xl_dict["mosaicity"] = mosaicity

        return xl_dict

    @classmethod
    def from_dict(cls, d):
        """Convert the dictionary to a crystal model

        Params:
            d The dictionary of parameters

        Returns:
            The crystal model

        """
        xl = MosaicCrystalKabsch2010(super(MosaicCrystalKabsch2010, cls).from_dict(d))

        # This parameter doesn't survive the Crystal copy constructor so has to be re-set.
        # Isoforms used for stills
        try:
            xl.identified_isoform = d["identified_isoform"]
        except KeyError:
            pass

        # Extract mosaicity, if present
        try:
            mosaicity = d["mosaicity"]
            xl.set_mosaicity(mosaicity)
        except KeyError:
            pass

        return xl


@boost.python.inject_into(MosaicCrystalSauter2014)
class _(object):
    def as_str(self, show_scan_varying=False):
        return "\n".join(
            (
                super(MosaicCrystalSauter2014, self).as_str(
                    show_scan_varying=show_scan_varying
                ),
                "    Half mosaic angle (degrees):  %.6f"
                % self.get_half_mosaicity_deg(),
                "    Domain size (Angstroms):  %.6f" % self.get_domain_size_ang(),
            )
        )

    def get_A_as_sqr(self):  # required for lunus
        return matrix.sqr(self.get_A())

    def get_A_inverse_as_sqr(self):
        return self.get_A_as_sqr().inverse()

    def to_dict(self):
        """Convert the crystal model to a dictionary

        Returns:
            A dictionary of the parameters

        """
        xl_dict = super(MosaicCrystalSauter2014, self).to_dict()

        # Get the mosaic parameters
        half_mosaicity = self.get_half_mosaicity_deg()
        xl_dict["ML_half_mosaicity_deg"] = half_mosaicity

        domain_size = self.get_domain_size_ang()
        xl_dict["ML_domain_size_ang"] = domain_size

        return xl_dict

    @classmethod
    def from_dict(cls, d):
        """Convert the dictionary to a crystal model

        Params:
            d The dictionary of parameters

        Returns:
            The crystal model

        """
        xl = MosaicCrystalSauter2014(super(MosaicCrystalSauter2014, cls).from_dict(d))

        # Parameters for maximum likelihood values
        try:
            xl.set_half_mosaicity_deg(d["ML_half_mosaicity_deg"])
        except KeyError:
            pass
        try:
            xl.set_domain_size_ang(d["ML_domain_size_ang"])
        except KeyError:
            pass

        # This parameter doesn't survive the Crystal copy constructor so have to be re-set
        # Isoforms used for stills
        try:
            xl.identified_isoform = d["identified_isoform"]
        except KeyError:
            pass

        return xl


@boost.python.inject_into(Experiment)
class _(object):
    def load_models(self, index=None):
        """ Load the models from the imageset """
        if index is None:
            index = 0
        self.beam = self.imageset.get_beam(index)
        self.detector = self.imageset.get_detector(index)
        self.goniometer = self.imageset.get_goniometer(index)
        self.scan = self.imageset.get_scan(index)


@boost.python.inject_into(ExperimentList)
class _(object):
    def __repr__(self):
        if len(self):
            return "ExperimentList([{}])".format(", ".join(repr(x) for x in self))
        else:
            return "ExperimentList()"

    def beams(self):
        """ Get a list of the unique beams (includes None). """
        return list(OrderedSet(e.beam for e in self))

    def detectors(self):
        """ Get a list of the unique detectors (includes None). """
        return list(OrderedSet(e.detector for e in self))

    def goniometers(self):
        """ Get a list of the unique goniometers (includes None). """
        return list(OrderedSet(e.goniometer for e in self))

    def scans(self):
        """ Get a list of the unique scans (includes None). """
        return list(OrderedSet(e.scan for e in self))

    def crystals(self):
        """ Get a list of the unique crystals (includes None). """
        return list(OrderedSet(e.crystal for e in self))

    def profiles(self):
        """ Get a list of the unique profile models (includes None). """
        return list(OrderedSet(e.profile for e in self))

    def scaling_models(self):
        """ Get a list of the unique scaling models (includes None). """
        return list(OrderedSet(e.scaling_model for e in self))

    def imagesets(self):
        """Get a list of the unique imagesets."""
        return list(OrderedSet([e.imageset for e in self if e.imageset is not None]))

    def all_stills(self):
        """Check if all the experiments are stills"""
        return all(exp.is_still() for exp in self)

    def all_sequences(self):
        """Check if all the experiments are from sequences"""
        return all(exp.is_sequence() for exp in self)

    def to_dict(self):
        """ Serialize the experiment list to dictionary. """

        # Check the experiment list is consistent
        assert self.is_consistent()

        # Table of names, source collections and experiment accessors
        # for member models that are serialized as an index table-lookup
        lookup_members = [
            ("beam", self.beams, lambda x: x.beam),
            ("detector", self.detectors, lambda x: x.detector),
            ("goniometer", self.goniometers, lambda x: x.goniometer),
            ("scan", self.scans, lambda x: x.scan),
            ("crystal", self.crystals, lambda x: x.crystal),
            ("profile", self.profiles, lambda x: x.profile),
            ("scaling_model", self.scaling_models, lambda x: x.scaling_model),
            ("imageset", self.imagesets, lambda x: x.imageset),
        ]
        # Generate index lookup tables for each output member collection instance
        index_lookup = {
            name: collections.OrderedDict(
                [
                    (model, i)
                    for i, model in enumerate(x for x in models() if x is not None)
                ]
            )
            for name, models, _ in lookup_members
        }

        # Create the output dictionary
        result = collections.OrderedDict()
        result["__id__"] = "ExperimentList"
        result["experiment"] = []

        # Add the experiments to the dictionary
        for e in self:
            obj = collections.OrderedDict()
            obj["__id__"] = "Experiment"
            obj["identifier"] = e.identifier

            # For each member model, look up the index
            for name, _, attr in lookup_members:
                model = attr(e)
                if model is not None:
                    obj[name] = index_lookup[name][model]

            result["experiment"].append(obj)

        def get_template(imset):
            if imset.reader().is_single_file_reader():
                return imset.reader().master_path()
            else:
                return imset.get_template()

        # Serialize all the imagesets
        result["imageset"] = []
        for imset in index_lookup["imageset"]:
            if isinstance(imset, ImageSequence):
                # FIXME_HACK
                template = get_template(imset)
                r = collections.OrderedDict(
                    [("__id__", "ImageSequence"), ("template", template)]
                )
                # elif isinstance(imset, MemImageSet):
                #   r = collections.OrderedDict([
                #     ('__id__', 'MemImageSet')])
                if imset.reader().is_single_file_reader():
                    r["single_file_indices"] = list(imset.indices())
            elif isinstance(imset, ImageSet):
                r = collections.OrderedDict(
                    [("__id__", "ImageSet"), ("images", imset.paths())]
                )
                if imset.reader().is_single_file_reader():
                    r["single_file_indices"] = list(imset.indices())
            elif isinstance(imset, ImageGrid):
                r = collections.OrderedDict(
                    [
                        ("__id__", "ImageGrid"),
                        ("images", imset.paths()),
                        ("grid_size", imset.get_grid_size()),
                    ]
                )
                if imset.reader().is_single_file_reader():
                    r["single_file_indices"] = list(imset.indices())
            else:
                raise TypeError(
                    "expected ImageSet or ImageSequence, got %s" % type(imset)
                )
            r["mask"] = imset.external_lookup.mask.filename
            r["gain"] = imset.external_lookup.gain.filename
            r["pedestal"] = imset.external_lookup.pedestal.filename
            r["dx"] = imset.external_lookup.dx.filename
            r["dy"] = imset.external_lookup.dy.filename
            r["params"] = imset.params()
            result["imageset"].append(r)

        # Extract all the ordered model dictionaries - is important these
        # preserve the same order as used in experiment serialization above
        for name, models in index_lookup.items():
            # Only fill out entries not handled above e.g. imageset
            if name not in result:
                result[name] = [x.to_dict() for x in models]

        # Return the dictionary
        return result

    def to_datablocks(self):
        """Return the experiment list as a datablock list.
        This assumes that the experiment contains 1 datablock."""
        # Datablock depends on model/__init__
        from dxtbx.datablock import DataBlockFactory

        # Convert the experiment list to dict
        obj = self.to_dict()
        # Convert the dictionary to a datablock dictionary
        obj["__id__"] = "DataBlock"
        for e in obj["experiment"]:
            iid = e["imageset"]
            imageset = obj["imageset"][iid]
            if "beam" in e:
                imageset["beam"] = e["beam"]
            if "detector" in e:
                imageset["detector"] = e["detector"]
            if "goniometer" in e:
                imageset["goniometer"] = e["goniometer"]
            if "scan" in e:
                imageset["scan"] = e["scan"]

        # Remove the experiments
        del obj["experiment"]

        # Create the datablock
        return DataBlockFactory.from_dict([obj])

    def nullify_all_single_file_reader_format_instances(self):
        """
        Parallel reading of HDF5 from the same handle is not allowed. Python
        multiprocessing is a bit messed up and used fork on linux so need to
        close and reopen file.

        """
        for experiment in self:
            if experiment.imageset.reader().is_single_file_reader():
                experiment.imageset.reader().nullify_format_instance()

    def as_json(self, filename=None, compact=False, split=False):
        """ Dump experiment list as json """
        # Get the dictionary and get the JSON string
        dictionary = self.to_dict()

        # Split into separate files
        if filename is not None and split:

            # Get lists of models by filename
            basepath = os.path.splitext(filename)[0]
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

        # Datablock depends on model/__init__
        from dxtbx.datablock import AutoEncoder

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
            if fname:
                with open(fname, "w") as outfile:
                    outfile.write(text)
            else:
                return text

    def as_pickle(self, filename=None, **kwargs):
        """ Dump experiment list as pickle. """
        # Get the pickle string
        text = pickle.dumps(self, protocol=pickle.HIGHEST_PROTOCOL)

        # Write the file
        if filename:
            with open(str(filename), "wb") as outfile:
                outfile.write(text)
        else:
            return text

    def as_file(self, filename, **kwargs):
        """ Dump experiment list as file. """
        ext = os.path.splitext(filename)[1]
        j_ext = [".json", ".expt"]
        p_ext = [".p", ".pkl", ".pickle"]
        if ext.lower() in j_ext:
            return self.as_json(filename, **kwargs)
        elif ext.lower() in p_ext:
            return self.as_pickle(filename, **kwargs)
        else:
            ext_str = "|".join(j_ext + p_ext)
            raise RuntimeError("expected extension {%s}, got %s" % (ext_str, ext))

    @staticmethod
    def from_file(filename, check_format=True):
        # type: (str, bool) -> ExperimentList
        """
        Load an ExperimentList from a serialized file.

        Args:
            filename: The filename to load an ExperimentList from
            check_format: If True, will attempt to verify image data type
        """
        # Inline to avoid recursive imports
        from .experiment_list import ExperimentListFactory

        return ExperimentListFactory.from_serialized_format(str(filename), check_format)


@boost.python.inject_into(Beam)
class _(object):
    def get_direction(self):
        warnings.warn(
            "Calling get_direction is deprecated. Please use "
            ".get_sample_to_source_direction() instead. "
            "See https://github.com/cctbx/dxtbx/issues/6",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_sample_to_source_direction()
