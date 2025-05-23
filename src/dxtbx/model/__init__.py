from __future__ import annotations

import copy
import json
import os
import sys

from ordered_set import OrderedSet

import boost_adaptbx.boost.python
import cctbx.crystal
import cctbx.sgtbx
import cctbx.uctbx
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.imageset import ImageGrid, ImageSequence, ImageSet
from dxtbx.model.beam import BeamFactory
from dxtbx.model.crystal import CrystalFactory
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model.profile import ProfileModelFactory
from dxtbx.model.scan import ScanFactory
from dxtbx.util import AutoEncoder, format_float_with_standard_uncertainty

try:
    from ..dxtbx_model_ext import (
        Beam,
        BeamBase,
        Crystal,
        CrystalBase,
        Detector,
        DetectorNode,
        Experiment,
        ExperimentList,
        ExperimentType,
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
        PolychromaticBeam,
        PxMmStrategy,
        Scan,
        ScanBase,
        SimplePxMmStrategy,
        Spectrum,
        VirtualPanel,
        VirtualPanelFrame,
        get_mod2pi_angles_in_range,
        get_range_of_mod2pi_angles,
        is_angle_in_range,
        parallax_correction,
        parallax_correction_inv,
    )
except ModuleNotFoundError:
    from dxtbx_model_ext import (  # type: ignore
        Beam,
        BeamBase,
        Crystal,
        CrystalBase,
        Detector,
        DetectorNode,
        Experiment,
        ExperimentList,
        ExperimentType,
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
        PolychromaticBeam,
        PxMmStrategy,
        Scan,
        ScanBase,
        SimplePxMmStrategy,
        Spectrum,
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
    "PolychromaticBeam",
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
    "Spectrum",
    "VirtualPanel",
    "VirtualPanelFrame",
    "get_mod2pi_angles_in_range",
    "get_range_of_mod2pi_angles",
    "is_angle_in_range",
    "parallax_correction",
    "parallax_correction_inv",
)


@boost_adaptbx.boost.python.inject_into(Detector)
class _detector:
    def iter_panels(self):
        """Iterate through just the panels depth-first."""
        for obj in self.iter_preorder():
            if obj.is_panel():
                yield obj

    def iter_preorder(self):
        """Iterate through the groups and panels depth-first."""
        stack = [self.hierarchy()]
        while stack:
            node = stack.pop()
            yield node
            if node.is_group():
                stack.extend(reversed(node))

    def iter_levelorder(self):
        """Iterate through the groups and panels breadth-first."""
        queue = [self.hierarchy()]
        while queue:
            node = queue.pop(0)
            yield node
            if node.is_group():
                queue.extend(node)


@boost_adaptbx.boost.python.inject_into(Crystal)
class _crystal:
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
            .mathematica_form(format="% 7.6f", one_row_per_line=True)
            .splitlines()
        )
        bmat = (
            matrix.sqr(self.get_B())
            .mathematica_form(format="% 7.6f", one_row_per_line=True)
            .splitlines()
        )
        amat = (
            (matrix.sqr(self.get_U()) * matrix.sqr(self.get_B()))
            .mathematica_form(format="% 7.6f", one_row_per_line=True)
            .splitlines()
        )

        msg = ["Crystal:"]

        if len(uc_sd) != 0:
            cell_str = [
                format_float_with_standard_uncertainty(v, e, minimum=1.0e-5)
                for (v, e) in zip(uc, uc_sd)
            ]
            msg.append("    Unit cell: " + ", ".join(cell_str))
        else:
            msg.append(
                "    Unit cell: " + "%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f" % uc
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
                        + "%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f" % uc
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

        uc = self.get_recalculated_unit_cell()
        if uc is not None:
            uc = uc.parameters()
            uc_sd = self.get_recalculated_cell_parameter_sd()
            if len(uc_sd) != 0:
                cell_str = [
                    format_float_with_standard_uncertainty(v, e, minimum=1.0e-5)
                    for (v, e) in zip(uc, uc_sd)
                ]
                msg.append("    Recalculated unit cell: " + ", ".join(cell_str))
            else:
                msg.append(
                    "    Recalculated unit cell: "
                    + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc
                )
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
        xl_dict = {
            "__id__": "crystal",
            "real_space_a": real_space_a,
            "real_space_b": real_space_b,
            "real_space_c": real_space_c,
            "space_group_hall_symbol": hall,
        }

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

        recalculated_unit_cell = self.get_recalculated_unit_cell()
        if recalculated_unit_cell is not None:
            xl_dict["recalculated_unit_cell"] = recalculated_unit_cell.parameters()
            xl_dict["recalculated_cell_parameter_sd"] = (
                self.get_recalculated_cell_parameter_sd()
            )
            xl_dict["recalculated_cell_volume_sd"] = (
                self.get_recalculated_cell_volume_sd()
            )

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

        recalculated_unit_cell = d.get("recalculated_unit_cell")
        if recalculated_unit_cell is not None:
            xl.set_recalculated_unit_cell(cctbx.uctbx.unit_cell(recalculated_unit_cell))

        recalculated_cell_parameter_sd = d.get("recalculated_cell_parameter_sd")
        if recalculated_cell_parameter_sd is not None:
            xl.set_recalculated_cell_parameter_sd(recalculated_cell_parameter_sd)

        recalculated_cell_volume_sd = d.get("recalculated_cell_volume_sd")
        if recalculated_cell_volume_sd is not None:
            xl.set_recalculated_cell_volume_sd(recalculated_cell_volume_sd)

        return xl


@boost_adaptbx.boost.python.inject_into(MosaicCrystalKabsch2010)
class _crystal_kabsch:
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


@boost_adaptbx.boost.python.inject_into(MosaicCrystalSauter2014)
class _crystalsauter:
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


@boost_adaptbx.boost.python.inject_into(Experiment)
class _experiment:
    def load_models(self, index=None):
        """Load the models from the imageset"""
        if index is None:
            index = 0
        self.beam = self.imageset.get_beam(index)
        self.detector = self.imageset.get_detector(index)
        self.goniometer = self.imageset.get_goniometer(index)
        self.scan = self.imageset.get_scan(index)


@boost_adaptbx.boost.python.inject_into(ExperimentList)
class _experimentlist:
    def __repr__(self):
        if len(self):
            return "ExperimentList([{}])".format(", ".join(repr(x) for x in self))
        else:
            return "ExperimentList()"

    def beams(self):
        """Get a list of the unique beams (includes None)."""
        return list(OrderedSet(e.beam for e in self))

    def detectors(self):
        """Get a list of the unique detectors (includes None)."""
        return list(OrderedSet(e.detector for e in self))

    def goniometers(self):
        """Get a list of the unique goniometers (includes None)."""
        return list(OrderedSet(e.goniometer for e in self))

    def scans(self):
        """Get a list of the unique scans (includes None)."""
        return list(OrderedSet(e.scan for e in self))

    def crystals(self):
        """Get a list of the unique crystals (includes None)."""
        return list(OrderedSet(e.crystal for e in self))

    def profiles(self):
        """Get a list of the unique profile models (includes None)."""
        return list(OrderedSet(e.profile for e in self))

    def scaling_models(self):
        """Get a list of the unique scaling models (includes None)."""
        return list(OrderedSet(e.scaling_model for e in self))

    def imagesets(self):
        """Get a list of the unique imagesets."""
        return list(OrderedSet([e.imageset for e in self if e.imageset is not None]))

    def all_stills(self):
        """Check if all the experiments are stills"""
        return all(exp.get_type() == ExperimentType.STILL for exp in self)

    def all_sequences(self):
        """Check if all the experiments are from sequences"""
        return self.all_rotations()

    def all_rotations(self):
        """Check if all the experiments are stills"""
        return all(exp.get_type() == ExperimentType.ROTATION for exp in self)

    def all_tof(self):
        """Check if all the experiments are time-of-flight"""
        return all(exp.get_type() == ExperimentType.TOF for exp in self)

    def all_laue(self):
        """Check if all the experiments are Laue experiments"""
        return all(exp.get_type() == ExperimentType.LAUE for exp in self)

    def all_same_type(self):
        """Check if all experiments are the same type"""
        if len(self) <= 1:
            return True
        expt_type = self[0].get_type()
        for i in range(1, len(self)):
            if self[i].get_type() != expt_type:
                return False
        return True

    def to_dict(self):
        """Serialize the experiment list to dictionary."""

        def abspath_or_none(filename):
            if filename is None or filename == "":
                return None
            return os.path.abspath(filename)

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
            name: {
                model: i for i, model in enumerate(x for x in models() if x is not None)
            }
            for name, models, _ in lookup_members
        }

        # Create the output dictionary
        result = {
            "__id__": "ExperimentList",
            "experiment": [],
        }

        # Add the experiments to the dictionary
        for e in self:
            obj = {
                "__id__": "Experiment",
                "identifier": e.identifier,
            }

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
                r = {
                    "__id__": "ImageSequence",
                    "template": template,
                }
                if imset.reader().is_single_file_reader():
                    r["single_file_indices"] = list(imset.indices())
            elif isinstance(imset, ImageSet):
                r = {
                    "__id__": "ImageSet",
                    "images": imset.paths(),
                }
                if imset.reader().is_single_file_reader():
                    r["single_file_indices"] = list(imset.indices())
            elif isinstance(imset, ImageGrid):
                r = {
                    "__id__": "ImageGrid",
                    "images": imset.paths(),
                    "grid_size": imset.get_grid_size(),
                }
                if imset.reader().is_single_file_reader():
                    r["single_file_indices"] = list(imset.indices())
            else:
                raise TypeError(
                    "expected ImageSet or ImageSequence, got %s" % type(imset)
                )
            r["mask"] = abspath_or_none(imset.external_lookup.mask.filename)
            r["gain"] = abspath_or_none(imset.external_lookup.gain.filename)
            r["pedestal"] = abspath_or_none(imset.external_lookup.pedestal.filename)
            r["dx"] = abspath_or_none(imset.external_lookup.dx.filename)
            r["dy"] = abspath_or_none(imset.external_lookup.dy.filename)
            r["params"] = imset.params()
            result["imageset"].append(r)

        # Extract all the ordered model dictionaries - is important these
        # preserve the same order as used in experiment serialization above
        for name, models in index_lookup.items():
            # Only fill out entries not handled above e.g. imageset
            if name not in result:
                result[name] = [x.to_dict() for x in models]

        return result

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
        """Dump experiment list as json"""
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
            edict = {
                "__id__": "ExperimentList",
                "experiment": dictionary["experiment"],
            }

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
            if fname:
                with open(fname, "w") as outfile:
                    outfile.write(text)
            else:
                return text

    def as_file(self, filename, **kwargs):
        """Dump experiment list as file."""
        ext = os.path.splitext(filename)[1]
        j_ext = [".json", ".expt"]
        if ext.lower() in j_ext:
            return self.as_json(filename, **kwargs)
        else:
            ext_str = "|".join(j_ext)
            raise RuntimeError(f"expected extension {{{ext_str}}}, got {ext}")

    @staticmethod
    def from_file(filename: str, check_format: bool = True) -> ExperimentList:
        """
        Load an ExperimentList from a serialized file.

        Args:
            filename: The filename to load an ExperimentList from
            check_format: If True, will attempt to verify image data type
        """
        # Inline to avoid recursive imports
        from .experiment_list import ExperimentListFactory

        return ExperimentListFactory.from_serialized_format(filename, check_format)

    def change_basis(self, change_of_basis_ops, in_place=False):
        """
        Apply change of basis operators to an ExperimentList

        Args:
            change_of_basis_ops: This can either be a single
                cctbx.sgtbx.change_of_basis_op to be applied to all experiments, or a
                list of operators, one per experiment.
            in_place (bool): Apply the change of basis operations in-place to the
                current ExperimentList. Default is to return a copy of the
                ExperimentList.

        Returns:
            The reindexed ExperimentList
        """
        if isinstance(change_of_basis_ops, cctbx.sgtbx.change_of_basis_op):
            change_of_basis_ops = [change_of_basis_ops] * len(self)
        assert len(change_of_basis_ops) == len(self), (
            "Number of change_of_basis_ops (%i) not equal to the number of experiments (%i)"
            % (len(change_of_basis_ops), len(self))
        )
        if in_place:
            return_expts = self
        else:
            return_expts = copy.deepcopy(self)
        for expt, cb_op in zip(return_expts, change_of_basis_ops):
            expt.crystal = expt.crystal.change_basis(cb_op)
        return return_expts

    @staticmethod
    def from_templates(templates, **kwargs):
        """Import an experiment list from templates"""
        # Import here to avoid cyclic dependencies
        from .experiment_list import ExperimentListFactory

        return ExperimentListFactory.from_templates(templates, **kwargs)
