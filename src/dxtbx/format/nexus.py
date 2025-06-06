from __future__ import annotations

import collections
import itertools
import math
import os
from collections.abc import Iterable
from typing import Union

import h5py
import hdf5plugin  # noqa: F401
import numpy

import cctbx.uctbx
from cctbx.eltbx import attenuation_coefficient
from scitbx import matrix
from scitbx.array_family import flex
from scitbx.matrix import col, sqr

import dxtbx.model
from dxtbx.model import (
    Beam,
    Crystal,
    Detector,
    Panel,
    ParallaxCorrectedPxMmStrategy,
    Scan,
    Spectrum,
)

try:
    try:
        from ..dxtbx_format_nexus_ext import (
            dataset_as_flex_double,
            dataset_as_flex_float,
            dataset_as_flex_int,
        )
    except ModuleNotFoundError:
        from dxtbx_format_nexus_ext import (  # type: ignore
            dataset_as_flex_double,
            dataset_as_flex_float,
            dataset_as_flex_int,
        )
except ImportError:
    # Workaround for psana build, which doesn't link HDF5 properly
    if "SIT_ROOT" not in os.environ:
        raise

NXNode = Union[h5py.File, h5py.Group]


def h5str(h5_value: str | numpy.bytes_ | bytes) -> str:
    """
    Convert a value returned an h5py attribute to str.

    h5py can return either a bytes-like (numpy.bytes_) or str object
    for attribute values depending on whether the value was written as
    fixed or variable length. This function collapses the two to str.
    """
    if isinstance(h5_value, (numpy.bytes_, bytes)):
        return h5_value.decode("utf-8")
    return h5_value


def dataset_as_flex(dataset, selection):
    if numpy.issubdtype(dataset.dtype, numpy.integer):
        return dataset_as_flex_int(dataset.id.id, selection)
    else:
        assert numpy.issubdtype(dataset.dtype, numpy.floating)
        double_types = [
            numpy.double,
            numpy.longdouble,
            numpy.float64,
        ]
        if hasattr(numpy, "float96"):
            double_types.append(numpy.float96)
        if hasattr(numpy, "float128"):
            double_types.append(numpy.float128)
        if dataset.dtype in [
            numpy.half,
            numpy.single,
            numpy.float64,
            numpy.float16,
            numpy.float32,
        ]:
            return dataset_as_flex_float(dataset.id.id, selection)
        elif dataset.dtype in double_types:
            return dataset_as_flex_double(dataset.id.id, selection)
        else:
            assert False, "unknown floating data type (%s)" % str(dataset.dtype)


class NXValidationError(RuntimeError):
    """A specific exception to record validation errors"""


def find_classes(node: NXNode, *nx_classes: str | None) -> tuple[list[h5py.Group], ...]:
    """
    Find instances of multiple NXclass types within the children of the current node.

    Args:
        node: The input h5py node (h5py.File or h5py.Group).
        nx_classes: Names of NXclass types to search for.  If None, search for children
            without an NXclass.

    Returns:
        A list of matching nodes for each of the specified NX_class types.
    """
    results: dict[str | None, list[h5py.Group]] = {
        nx_class: [] for nx_class in nx_classes
    }

    values: Iterable[h5py.Group] = filter(None, node.values())
    for v in values:
        class_name = h5str(v.attrs.get("NX_class"))
        if class_name in nx_classes:
            results[class_name].append(v)

    return tuple(results.values())


def find_class(node: NXNode, nx_class: str | None) -> list[h5py.Group]:
    """
    Find instances of a single NXclass type within the children of the current node.

    This is a convenience function, equivalent to calling find_classes with a single
    NXclass type name argument and returning the list of matches.

    Args:
        node: The input h5py node (h5py.File or h5py.Group).
        nx_class: Names of NXclass type to search for.  If None, search for children
            without an NXclass.

    Returns:
        The list of matching nodes for the specified NXclass type.
    """
    return find_classes(node, nx_class)[0]


def find_entries(nx_file: h5py.File) -> list[h5py.Group]:
    """
    Find NXmx entries.

    A valid NeXus file should contain one or more NXentry groups at the top level. Only
    return those NXentry groups that are specified as following the NXmx definition.

    https://manual.nexusformat.org/classes/base_classes/NXentry.html#nxentry

    Args:
        nx_file (h5py.File): The input h5py.File instance.

    Returns:
        list: The list of found NXentry groups.

    """
    return [
        group
        for group in find_class(nx_file, "NXentry")
        if "definition" in group and h5str(group["definition"][()]) == "NXmx"
    ]


def convert_units(value, input_units, output_units):
    """
    Hacky utility function to convert units
    """
    if isinstance(input_units, bytes):
        input_units = input_units.decode("latin-1")
    converters = {
        "m": {
            "mm": lambda x: x * 1e3,
            "microns": lambda x: x * 1e6,
            "nm": lambda x: x * 1e9,
        },
        "mm": {
            "m": lambda x: x * 1e-3,
            "microns": lambda x: x * 1e3,
            "nm": lambda x: x * 1e6,
        },
        "microns": {
            "m": lambda x: x * 1e-6,
            "mm": lambda x: x * 1e-3,
            "nm": lambda x: x * 1e3,
        },
        "nm": {
            "m": lambda x: x * 1e-9,
            "mm": lambda x: x * 1e-6,
            "microns": lambda x: x * 1e-3,
            "angstroms": lambda x: x * 10,
        },
        "angstroms": {"angstrom": lambda x: x},
    }
    if input_units == output_units:
        return value
    try:
        return converters[input_units][output_units](value)
    except Exception:
        raise RuntimeError(f"Can't convert units {input_units!r} to {output_units!r}")


def visit_dependencies(nx_file, item, visitor=None):
    """
    Walk the dependency chain and call a visitor function
    """
    dependency_chain = set()
    if os.path.basename(item) == "depends_on":
        depends_on = nx_file[item][()]
    else:
        depends_on = h5str(nx_file[item].attrs["depends_on"])

    while not depends_on == ".":
        if visitor:
            visitor(nx_file, depends_on)
        if depends_on in dependency_chain:
            raise RuntimeError("'%s' is a circular dependency" % depends_on)
        try:
            _ = nx_file[depends_on]
        except Exception:
            raise RuntimeError("'%s' is missing from nx_file" % depends_on)
        dependency_chain.add(depends_on)
        try:
            depends_on = h5str(nx_file[depends_on].attrs["depends_on"])
        except Exception:
            raise RuntimeError("'%s' contains no depends_on attribute" % depends_on)


def construct_vector(nx_file, item, vector=None):
    """
    Walk the dependency chain and create the absolute vector
    """

    class TransformVisitor:
        def __init__(self, vector):
            self.vector = matrix.col(vector)

        def visit(self, nx_file, depends_on):
            item = nx_file[depends_on]
            value = item[()]
            units = h5str(item.attrs["units"])
            ttype = h5str(item.attrs["transformation_type"])
            vector = matrix.col(item.attrs["vector"])
            if ttype == "translation":
                value = convert_units(value, units, "mm")
                if hasattr(value, "__iter__") and len(value) == 1:
                    value = value[0]
                self.vector = vector * value + self.vector
            elif ttype == "rotation":
                if hasattr(value, "__iter__") and len(value):
                    value = value[0]
                if units == "rad":
                    deg = False
                elif units == "deg":
                    deg = True
                else:
                    raise RuntimeError("Invalid units: %s" % units)
                self.vector.rotate_around_origin(axis=vector, angle=value, deg=deg)
            else:
                raise RuntimeError("Unknown transformation_type: %s" % ttype)

    if vector is None:
        value = nx_file[item][()]
        units = h5str(nx_file[item].attrs["units"])
        ttype = h5str(nx_file[item].attrs["transformation_type"])
        vector = nx_file[item].attrs["vector"]
        if "offset" in nx_file[item].attrs:
            offset = nx_file[item].attrs["offset"]
            offset = convert_units(offset, units, "mm")
        else:
            offset = vector * 0.0
        if ttype == "translation":
            value = convert_units(value, units, "mm")
            try:
                vector = vector * value
            except ValueError:
                vector = vector * value[0]
            vector += offset

    visitor = TransformVisitor(vector)
    visit_dependencies(nx_file, item, visitor.visit)
    return visitor.vector


def construct_axes(nx_file, item, vector=None):
    """
    Walk the dependency chain and create the absolute vector
    """

    class Visitor:
        def __init__(self):
            self._axes = flex.vec3_double()
            self._angles = flex.double()
            self._axis_names = flex.std_string()
            self._is_scan_axis = flex.bool()

        def visit(self, nx_file, depends_on):
            item = nx_file[depends_on]
            value = item[()]
            units = h5str(item.attrs["units"])
            ttype = h5str(item.attrs["transformation_type"])
            vector = [float(v) for v in item.attrs["vector"]]
            if ttype == "translation":
                return
            elif ttype == "rotation":
                if hasattr(value, "__iter__") and len(value):
                    value = value[0]
                if units == "rad":
                    value *= 180 / math.pi
                elif units not in ["deg", "degree", "degrees"]:
                    raise RuntimeError("Invalid units: %s" % units)

                # is the axis moving? Check the values for this axis
                v = item[...]
                if v.min() < v.max():
                    is_scan_axis = True
                else:
                    is_scan_axis = False

                # Is different coordinate system called mcstas
                # Rotate 180 about up if memory serves
                axis_name = item.name.split("/")[-1]
                self._axes.append(vector)
                self._angles.append(float(value))
                self._axis_names.append(str(axis_name))
                self._is_scan_axis.append(is_scan_axis)

            else:
                raise RuntimeError("Unknown transformation_type: %s" % ttype)

        def result(self):
            if self._is_scan_axis.count(True) == 0:
                # XXX not sure how best to handle this, but probably a still so no scan axis
                scan_axis = 0
            else:
                assert self._is_scan_axis.count(True) == 1, (
                    "Only one axis can be a scan axis: %s" % list(self._is_scan_axis)
                )
                scan_axis = flex.first_index(self._is_scan_axis, True)

            # Rotate 180 about up from McStas coordinate system
            cb_op = (-1, 0, 0, 0, 1, 0, 0, 0, -1)
            return cb_op * self._axes, self._angles, self._axis_names, scan_axis

    if vector is None:
        value = nx_file[item][()]
        units = h5str(nx_file[item].attrs["units"])
        ttype = h5str(nx_file[item].attrs["transformation_type"])
        vector = nx_file[item].attrs["vector"]
        if "offset" in nx_file[item].attrs:
            offset = nx_file[item].attrs["offset"]
        else:
            offset = vector * 0.0
        if ttype == "translation":
            value = convert_units(value, units, "mm")
            try:
                vector = vector * value
            except ValueError:
                vector = vector * value[0]
            vector += offset

    visitor = Visitor()
    visitor.visit(nx_file, item)

    visit_dependencies(nx_file, item, visitor.visit)

    return visitor.result()


class NXdetector_module:
    """
    A class to hold a handle to NXdetector_module
    """

    def __init__(self, handle):
        self.handle = handle


class NXdetector_group:
    """
    A class to hold a handle to NXdetector_group
    """

    def __init__(self, handle):
        self.handle = handle


class NXdetector:
    """
    A class to handle a handle to NXdetector
    """

    def __init__(self, handle):
        self.handle = handle

        # Find the NXdetector_modules
        self.modules = [
            NXdetector_module(entry)
            for entry in find_class(self.handle, "NXdetector_module")
        ]

        # Check we've got some stuff
        if not self.modules:
            raise NXValidationError("No NXdetector_module in %s" % self.handle.name)


class NXinstrument:
    """
    A class to hold a handle to NXinstrument
    """

    def __init__(self, handle):
        self.handle = handle

        # Find the NXdetector, any detector groups and the NXbeam
        detectors, detector_groups, beams = find_classes(
            self.handle, "NXdetector", "NXdetector_group", "NXbeam"
        )

        # Check we've got stuff
        if not detectors:
            raise NXValidationError("No NXdetector in %s" % self.handle.name)

        self.detectors = [NXdetector(detector) for detector in detectors]
        self.detector_groups = [NXdetector_group(group) for group in detector_groups]
        self.beams = [NXbeam(beam) for beam in beams]


class NXbeam:
    """
    A class to hold a handle to NXbeam
    """

    def __init__(self, handle):
        self.handle = handle


class NXsample:
    """
    A class to hold a handle to NXsample
    """

    def __init__(self, handle):
        self.handle = handle

        # Find the NXbeam
        self.beams = [NXbeam(beam) for beam in find_class(self.handle, "NXbeam")]


class NXdata:
    """
    A class to hold a handle to NXdata
    """

    def __init__(self, handle):
        self.handle = handle


class NXmxEntry:
    """
    A class to hold a handle to NXmx entries
    """

    def __init__(self, handle):
        self.handle = handle

        # Find the NXinstrument, NXsample & NXdata
        instruments, samples, data_sets = find_classes(
            self.handle, "NXinstrument", "NXsample", "NXdata"
        )

        # Check we've got some stuff
        if not instruments:
            raise NXValidationError("No NXinstrument in %s" % self.handle.name)
        if not samples:
            raise NXValidationError("No NXsample in %s" % self.handle.name)
        if not data_sets:
            raise NXValidationError("No NXdata in %s" % self.handle.name)

        self.instruments = [NXinstrument(instrument) for instrument in instruments]
        self.samples = [NXsample(sample) for sample in samples]
        self.data = [NXdata(data) for data in data_sets]


class NXmxReader:
    """
    A hacky class to read an NXmx file
    """

    def __init__(self, filename=None, handle=None):
        # Get the file handle
        if filename is not None:
            handle = h5py.File(filename, "r", swmr=True)

        # Find the NXmx entries
        self.entries = []
        for entry in find_entries(handle):
            self.entries.append(NXmxEntry(entry))

        # Check we've got some stuff
        if not self.entries:
            raise RuntimeError(
                "Error reading NXmxfile %r. No NXmx entries in file" % filename
            )

    def print_description(self):
        """
        Print a description of the NXmx file

        """
        print(" > Found %d NXmx entries" % len(self.entries))
        for entry in self.entries:
            handle = entry.handle
            instruments = entry.instruments
            samples = entry.samples
            print("  > %s" % handle.name)
            beams = []
            for instrument in instruments:
                handle = instrument.handle
                beams += instrument.beams
                detectors = instrument.detectors
                print("   > %s" % handle.name)
                for detector in detectors:
                    handle = detector.handle
                    modules = detector.modules
                    print("    > %s" % handle.name)
                    for module in modules:
                        handle = module.handle
                        print("     > %s" % handle.name)
            for sample in samples:
                handle = sample.handle
                beams += sample.beams
                print("   > %s" % handle.name)
            for beam in beams:
                handle = beam.handle
                print("    > %s" % handle.name)


def is_nexus_file(filename):
    """
    A hacky function to check if this is a nexus file
    """
    with h5py.File(filename, "r") as handle:
        # Find the NXmx entries
        return bool(find_entries(handle))


class BeamFactory:
    """
    A class to create a beam model from NXmx stuff
    """

    def __init__(self, obj):
        self.obj = obj
        self.model = None
        self.index = None
        self.spectrum = None

    def read_models(self, index=None):
        self.load_model(index)
        return self.model, self.spectrum

    def load_model(self, index=None):
        # Cached model
        if self.model is not None and index == self.index:
            return self.model

        # Get the items from the NXbeam class
        primary_key = "incident_wavelength"
        wavelength = self.obj.handle[primary_key]
        spectrum_wavelengths = wavelength
        spectrum_weights = self.obj.handle.get(primary_key + "_weights")
        if spectrum_weights is None:
            # Handle deprecation: https://github.com/nexusformat/definitions/issues/837
            spectrum_weights = self.obj.handle.get(primary_key + "_weight")

        # If the wavelength array does not represent spectra, look for spectra
        # in the variant chain
        variant_test = wavelength
        has_variant_spectra = False
        while spectrum_weights is None:
            if "variant" in variant_test.attrs:
                variant_key = variant_test.attrs["variant"]
                variant_wavelengths = self.obj.handle[variant_key]
                variant_weights = self.obj.handle.get(variant_key + "_weights")
                if variant_weights is None:
                    # Handle deprecation: https://github.com/nexusformat/definitions/issues/837
                    variant_weights = self.obj.handle.get(variant_key + "_weight")
                if variant_weights is None:
                    variant_test = variant_wavelengths  # Keep looking
                else:
                    # Found spectra
                    spectrum_wavelengths = variant_wavelengths
                    spectrum_weights = variant_weights  # cause while loop to end
                    has_variant_spectra = True
            else:
                break

        if index is None:
            index = 0
        self.index = index

        def get_wavelength(wavelength):
            if wavelength.shape in ((), (1,)):
                wavelength_value = wavelength[()]
            else:
                wavelength_value = wavelength[index]
            wavelength_units = wavelength.attrs["units"]
            wavelength_value = float(
                convert_units(wavelength_value, wavelength_units, "angstrom")
            )
            return wavelength_value

        if spectrum_weights is None:
            # Construct the beam model
            wavelength_value = get_wavelength(wavelength)
            self.model = Beam(direction=(0, 0, 1), wavelength=wavelength_value)
        else:
            self.model = Beam()
            self.model.set_direction((0, 0, 1))

            wavelength_units = spectrum_wavelengths.attrs["units"]

            if len(spectrum_wavelengths.shape) > 1:
                spectrum_wavelengths = spectrum_wavelengths[index]
            else:
                spectrum_wavelengths = spectrum_wavelengths[()]
            if len(spectrum_weights.shape) > 1:
                spectrum_weights = spectrum_weights[index]
            else:
                spectrum_weights = spectrum_weights[()]

            spectrum_wavelengths = convert_units(
                spectrum_wavelengths, wavelength_units, "angstrom"
            )
            spectrum_energies = cctbx.factor_ev_angstrom / spectrum_wavelengths
            self.spectrum = Spectrum(spectrum_energies, spectrum_weights)

            if has_variant_spectra:
                wavelength_value = get_wavelength(wavelength)
                self.model.set_wavelength(wavelength_value)
            else:
                self.model.set_wavelength(self.spectrum.get_weighted_wavelength())
        return self.model


def get_change_of_basis(transformation, include_setting=True):
    """
    Get the 4x4 homogenous coordinate matrix for a given NXtransformation.
    """
    # Change of basis to convert from NeXus to IUCr/ImageCIF convention
    n2i_cob = sqr((-1, 0, 0, 0, 1, 0, 0, 0, -1))

    axis_type = h5str(transformation.attrs["transformation_type"])

    vector = n2i_cob * col(transformation.attrs["vector"]).normalize()
    setting = transformation[0] if include_setting else 0
    units = h5str(transformation.attrs["units"])

    if "offset" in transformation.attrs:
        offset = n2i_cob * col(transformation.attrs["offset"])
        if "offset_units" in transformation.attrs:
            offset_units = transformation.attrs["offset_units"]
        else:
            offset_units = units
        offset = convert_units(offset, offset_units, "mm")
    else:
        offset = col((0, 0, 0))

    # 4x4 change of basis matrix (homogeneous coordinates)
    cob = None

    if axis_type == "rotation":
        if units == "rad":
            deg = False
        elif units in ["deg", "degree", "degrees"]:
            deg = True
        else:
            raise RuntimeError("Invalid units: %s" % units)
        r3 = vector.axis_and_angle_as_r3_rotation_matrix(setting, deg=deg)
        cob = sqr(
            (
                r3[0],
                r3[1],
                r3[2],
                offset[0],
                r3[3],
                r3[4],
                r3[5],
                offset[1],
                r3[6],
                r3[7],
                r3[8],
                offset[2],
                0,
                0,
                0,
                1,
            )
        )
    elif axis_type == "translation":
        setting = convert_units(setting, units, "mm")
        translation = offset + (vector * setting)
        cob = sqr(
            (
                1,
                0,
                0,
                translation[0],
                0,
                1,
                0,
                translation[1],
                0,
                0,
                1,
                translation[2],
                0,
                0,
                0,
                1,
            )
        )
    else:
        raise ValueError("Unrecognized transformation type: %s" % axis_type)

    return cob


def get_depends_on_chain_using_equipment_components(transformation):
    """
    Given an NXtransformation, find the list of dependencies that transformation
    has.  If there are multiple dependencies in a single 'level', as indicated
    by grouping them using equipment_component, then the dependency chain will
    skip the intermediate dependencies, listing only the first at each level.
    """
    chain = [transformation]
    current = transformation

    while True:
        parent_id = h5str(current.attrs["depends_on"])

        if parent_id == ".":
            return chain
        parent = current.parent[parent_id]

        if "equipment_component" in current.attrs:
            eq_comp = current.attrs["equipment_component"]
            parent_eq_comp = parent.attrs["equipment_component"]
            if eq_comp == parent_eq_comp:
                current = parent
                continue
        chain.append(parent)
        current = parent


def get_cumulative_change_of_basis(transformation):
    """Get the 4x4 homogenous coordinate matrix for a given NXtransformation,
    combining it with the change of basis matrices of parent transformations
    with the same equipment component as the given transformation.
    Returns (parent, change of basis matrix), where parent is None if the
    transformation's depends_on is ".".  Parent is the transformation that the
    top level transformation in this chain of transformations depends on.
    """

    cob = get_change_of_basis(transformation)

    parent_id = h5str(transformation.attrs["depends_on"])

    if parent_id == ".":
        return None, cob
    parent = transformation.parent[parent_id]

    if "equipment_component" in transformation.attrs:
        eq_comp = transformation.attrs["equipment_component"]
        parent_eq_comp = parent.attrs["equipment_component"]
        if eq_comp == parent_eq_comp:
            non_matching_parent, parent_cob = get_cumulative_change_of_basis(parent)
            return non_matching_parent, parent_cob * cob

    return parent, cob


class DetectorFactoryFromGroup:
    """
    A class to create a detector model from a NXdetector_group
    """

    def __init__(self, instrument, beam, idx=None):
        assert len(instrument.detector_groups) == 1, "Multiple detectors not supported"

        nx_group = instrument.detector_groups[0].handle
        group_names = nx_group["group_names"]
        group_parent = nx_group["group_parent"]

        # Verify the NXdetector objects specified by the detector group are present
        expected_detectors = []
        root_name = None
        for i, parent_id in enumerate(group_parent):
            assert parent_id in [
                -1,
                1,
            ], (
                "Hierarchy of detectors not supported. Hierarchy of module components within detector elements is supported"
            )

            if parent_id == -1:
                assert root_name is None, "Multiple roots not supported"
                root_name = group_names[i]
            else:
                expected_detectors.append(group_names[i].astype(str))

        assert root_name is not None, "Detector root not found"
        assert sorted(
            os.path.basename(d.handle.name) for d in instrument.detectors
        ) == sorted(expected_detectors), (
            "Mismatch between detector group names and detectors available"
        )

        root = None

        def set_frame(pg, transformation):
            """Function to set the local frame of a panel/panel group given an
            NXtransformation"""
            parent, cob = get_cumulative_change_of_basis(transformation)
            # set up the dxtbx d matrix.  Note use of homogenous coordinates.
            origin = matrix.col((cob * matrix.col((0, 0, 0, 1)))[0:3])
            fast = matrix.col((cob * matrix.col((1, 0, 0, 1)))[0:3]) - origin
            slow = matrix.col((cob * matrix.col((0, 1, 0, 1)))[0:3]) - origin

            pg.set_local_frame(fast.elems, slow.elems, origin.elems)

            name = str(os.path.basename(transformation.name))
            pg.set_name(name)

        self.model = Detector()
        for nx_detector in instrument.detectors:
            # Get the detector type
            if "type" in nx_detector.handle:
                detector_type = str(nx_detector.handle["type"][()])
            else:
                detector_type = "unknown"

            # get the set of leaf modules (handles nested modules)
            modules = [
                module
                for module in nx_detector.modules
                if not find_class(module.handle, "NXdetector_module")
            ]

            # depends_on field for a detector will have NxM entries in it, where
            # N = number of images and M = number of detector modules

            for module_number, nx_detector_module in enumerate(modules):
                # Get the depends_on starting point for this module and for this image
                panel_name = str(os.path.basename(nx_detector_module.handle.name))
                # image size stored slow to fast but dxtbx needs fast to slow
                image_size = tuple(
                    reversed(
                        list(map(int, nx_detector_module.handle["data_size"][-2:]))
                    )
                )

                # Get the trusted range of pixel values - if missing use
                # full range of signed 32 bit int
                try:
                    min_trusted_value = float(nx_detector.handle["underload_value"][()])
                except KeyError:
                    min_trusted_value = -0x7FFFFFFF

                try:
                    max_trusted_value = float(
                        nx_detector.handle["saturation_value"][()]
                    )
                except KeyError:
                    max_trusted_value = 0x7FFFFFFF

                trusted_range = min_trusted_value, max_trusted_value

                fast_pixel_direction_handle = nx_detector_module.handle[
                    "fast_pixel_direction"
                ]
                slow_pixel_direction_handle = nx_detector_module.handle[
                    "slow_pixel_direction"
                ]
                assert (
                    fast_pixel_direction_handle.attrs["depends_on"]
                    == slow_pixel_direction_handle.attrs["depends_on"]
                )
                depends_on = fast_pixel_direction_handle
                fast_pixel_direction_value = convert_units(
                    fast_pixel_direction_handle[0],
                    fast_pixel_direction_handle.attrs["units"],
                    "mm",
                )
                slow_pixel_direction_value = convert_units(
                    slow_pixel_direction_handle[0],
                    slow_pixel_direction_handle.attrs["units"],
                    "mm",
                )
                pixel_size = (
                    float(fast_pixel_direction_value),
                    float(slow_pixel_direction_value),
                )

                # Set up the hierarchical detector by iteraing through the dependencies,
                # starting at the root
                chain = get_depends_on_chain_using_equipment_components(depends_on)
                chain.pop(0)  # Do not want fast/slow pixel directions here
                pg = None
                for transform in reversed(chain):
                    name = str(os.path.basename(transform.name))
                    if (
                        pg is None
                    ):  # The first transform will be the root of the hierarchy
                        if root is None:
                            root = self.model.hierarchy()
                            set_frame(root, transform)
                        else:
                            assert root.get_name() == name, "Found multiple roots"
                        pg = root
                        continue

                    # subsequent times through the loop, pg will be the parent of the
                    # current transform
                    pg_names = [child.get_name() for child in pg]
                    if name in pg_names:
                        pg = pg[pg_names.index(name)]
                    else:
                        pg = pg.add_group()
                        set_frame(pg, transform)

                # pg is now this panel's parent
                p = pg.add_panel()
                fast = fast_pixel_direction_handle.attrs["vector"]
                fast = matrix.col([-fast[0], fast[1], -fast[2]])
                slow = slow_pixel_direction_handle.attrs["vector"]
                slow = matrix.col([-slow[0], slow[1], -slow[2]])
                cob = get_change_of_basis(
                    fast_pixel_direction_handle, include_setting=False
                ) * get_change_of_basis(
                    slow_pixel_direction_handle, include_setting=False
                )
                origin = matrix.col((cob * matrix.col((0, 0, 0, 1)))[0:3])

                p.set_local_frame(fast.elems, slow.elems, origin.elems)

                p.set_name(panel_name)
                p.set_pixel_size(pixel_size)
                p.set_image_size(image_size)
                p.set_type(detector_type)
                p.set_trusted_range(trusted_range)

                if "sensor_thickness" in nx_detector.handle:
                    # Get the detector thickness
                    thickness = nx_detector.handle["sensor_thickness"]
                    thickness_value = float(thickness[()])
                    thickness_units = thickness.attrs["units"]
                    thickness_value = float(
                        convert_units(thickness_value, thickness_units, "mm")
                    )
                    p.set_thickness(thickness_value)

                # Get the detector material
                if "sensor_material" in nx_detector.handle:
                    value = h5str(nx_detector.handle["sensor_material"][()])
                    material = {
                        "Si": "Si",
                        "Silicon": "Si",
                        "Sillicon": "Si",
                        "CdTe": "CdTe",
                        "GaAs": "GaAs",
                    }.get(value)
                    if not material:
                        raise RuntimeError("Unknown material: %s" % value)
                    p.set_material(material)

                    # Compute the attenuation coefficient.
                    # This will fail for undefined composite materials
                    # mu_at_angstrom returns cm^-1, but need mu in mm^-1
                    table = attenuation_coefficient.get_table(material)
                    wavelength = beam.get_wavelength()
                    mu = table.mu_at_angstrom(wavelength) / 10.0
                    p.set_mu(mu)

                if (
                    "sensor_thickness" in nx_detector.handle
                    and "sensor_material" in nx_detector.handle
                ):
                    p.set_px_mm_strategy(
                        ParallaxCorrectedPxMmStrategy(mu, thickness_value)
                    )


class DetectorFactory:
    """
    A class to create a detector model from NXmx stuff
    """

    def __init__(self, obj, beam, shape=None):
        # Get the handles
        nx_file = obj.handle.file
        nx_detector = obj.handle
        nx_module = obj.modules[0].handle

        # Get the detector name and type
        if "type" in nx_detector:
            detector_type = str(nx_detector["type"][()])
        else:
            detector_type = "unknown"
        detector_name = str(nx_detector.name)

        try:
            # underload_value: The lowest value at which pixels for this detector would be reasonably be measured.
            # https://manual.nexusformat.org/classes/applications/NXmx.html#nxmx-entry-instrument-detector-underload-value-field
            min_trusted_value = float(nx_detector["underload_value"][()])
        except KeyError:
            min_trusted_value = -0x7FFFFFFF

        try:
            # saturation_value: The value at which the detector goes into saturation.
            # Data above this value is known to be invalid.
            # https://manual.nexusformat.org/classes/applications/NXmx.html#nxmx-entry-instrument-detector-saturation-value-field
            max_trusted_value = float(nx_detector["saturation_value"][()])
        except KeyError:
            max_trusted_value = 0x7FFFFFFF

        trusted_range = min_trusted_value, max_trusted_value

        # Get the detector thickness
        thickness = nx_detector["sensor_thickness"]
        thickness_value = float(thickness[()])
        thickness_units = thickness.attrs["units"]
        thickness_value = float(convert_units(thickness_value, thickness_units, "mm"))

        # Get the detector material
        material = {
            "Si": "Si",
            "Silicon": "Si",
            "Sillicon": "Si",
            "CdTe": "CdTe",
            "GaAs": "GaAs",
        }.get(h5str(nx_detector["sensor_material"][()]))
        if not material:
            raise RuntimeError(
                "Unknown material: %s" % nx_detector["sensor_material"][()]
            )

        try:
            x_pixel = nx_detector["x_pixel_size"][()] * 1000.0
            y_pixel = nx_detector["y_pixel_size"][()] * 1000.0

            legacy_beam_x = float(x_pixel * nx_detector["beam_center_x"][()])
            legacy_beam_y = float(y_pixel * nx_detector["beam_center_y"][()])
        except KeyError:
            legacy_beam_x = 0
            legacy_beam_y = 0

        # Get the fast pixel size and vector
        fast_pixel_direction = nx_module["fast_pixel_direction"]
        fast_pixel_direction_value = float(fast_pixel_direction[()])
        fast_pixel_direction_units = fast_pixel_direction.attrs["units"]
        fast_pixel_direction_vector = fast_pixel_direction.attrs["vector"]
        fast_pixel_direction_value = convert_units(
            fast_pixel_direction_value, fast_pixel_direction_units, "mm"
        )
        fast_axis = matrix.col(fast_pixel_direction_vector).normalize()

        # Get the slow pixel size and vector
        slow_pixel_direction = nx_module["slow_pixel_direction"]
        slow_pixel_direction_value = float(slow_pixel_direction[()])
        slow_pixel_direction_units = slow_pixel_direction.attrs["units"]
        slow_pixel_direction_vector = slow_pixel_direction.attrs["vector"]
        slow_pixel_direction_value = convert_units(
            slow_pixel_direction_value, slow_pixel_direction_units, "mm"
        )
        slow_axis = matrix.col(slow_pixel_direction_vector).normalize()

        # Get the origin vector - working around if absent
        module_offset = nx_module["module_offset"]
        origin = construct_vector(nx_file, module_offset.name)
        if origin.elems[0] == 0.0 and origin.elems[1] == 0:
            origin = -(origin + legacy_beam_x * fast_axis + legacy_beam_y * slow_axis)

        # Change of basis to convert from NeXus to IUCr/ImageCIF convention
        cob = matrix.sqr((-1, 0, 0, 0, 1, 0, 0, 0, -1))
        origin = cob * matrix.col(origin)
        fast_axis = cob * fast_axis
        slow_axis = cob * slow_axis

        # Ensure that fast and slow axis are orthogonal
        normal = fast_axis.cross(slow_axis)
        slow_axis = -fast_axis.cross(normal)

        # Compute the attenuation coefficient.
        # This will fail for undefined composite materials
        # mu_at_angstrom returns cm^-1, but need mu in mm^-1
        table = attenuation_coefficient.get_table(material)
        wavelength = beam.get_wavelength()
        mu = float(table.mu_at_angstrom(wavelength)) / 10.0

        # Construct the detector model
        pixel_size = (fast_pixel_direction_value, slow_pixel_direction_value)

        # image size stored slow to fast but dxtbx needs fast to slow - assume
        # that the shapes being taken as input are slow -> fast
        if shape:
            image_size = tuple(reversed(shape[-2:]))
        else:
            image_size = tuple(int(x) for x in reversed(nx_module["data_size"][-2:]))

        self.model = Detector()
        self.model.add_panel(
            Panel(
                detector_type,
                detector_name,
                tuple(float(x) for x in fast_axis),
                tuple(float(x) for x in slow_axis),
                tuple(float(x) for x in origin),
                pixel_size,
                image_size,
                trusted_range,
                thickness_value,
                material,
                mu,
            )
        )

        # Set the parallax correction
        for panel in self.model:
            panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, thickness_value))
            panel.set_type("SENSOR_PAD")


class GoniometerFactory:
    """
    A class to create a goniometer model from NXmx stuff
    """

    def __init__(self, obj):
        if h5str(obj.handle["depends_on"][()]) == ".":
            self.model = None
        else:
            axes, angles, axis_names, scan_axis = construct_axes(
                obj.handle.file, obj.handle.file[obj.handle["depends_on"][()]].name
            )

            if len(axes) == 1:
                self.model = dxtbx.model.GoniometerFactory.make_goniometer(
                    axes[0], (1, 0, 0, 0, 1, 0, 0, 0, 1)
                )
            else:
                self.model = dxtbx.model.GoniometerFactory.make_multi_axis_goniometer(
                    axes, angles, axis_names, scan_axis
                )


def find_goniometer_rotation(obj):
    if h5str(obj.handle["depends_on"][()]) == ".":
        return
    thing = obj.handle.file[obj.handle["depends_on"][()]]
    tree = get_depends_on_chain_using_equipment_components(thing)
    for t in tree:
        o = obj.handle.file[t.name]
        if h5str(o.attrs["transformation_type"]) == "rotation":
            # if this is changing, assume is scan axis
            v = o[...]
            if v.min() < v.max():
                return o
    raise ValueError("no rotation found")


def find_scanning_axis(obj):
    if h5str(obj.handle["depends_on"][()]) == ".":
        return
    thing = obj.handle.file[obj.handle["depends_on"][()]]
    tree = get_depends_on_chain_using_equipment_components(thing)
    for t in tree:
        o = obj.handle.file[t.name]
        if o[()].size > 1:
            return o


def generate_scan_model(obj, detector_obj):
    """
    Create a scan model from NXmx stuff.
    """
    if h5str(obj.handle["depends_on"][()]) == ".":
        return

    # Get the image and oscillation range - need to search for rotations
    # in dependency tree - if not, find translations or just the thing
    # the sample depends on
    try:
        scan_axis = find_goniometer_rotation(obj)
    except ValueError:
        scan_axis = find_scanning_axis(obj)

    if scan_axis is None:
        scan_axis = obj.handle.file[obj.handle["depends_on"][()]]

    num_images = len(scan_axis)
    image_range = (1, num_images)

    rotn = h5str(scan_axis.attrs["transformation_type"]) == "rotation"

    if num_images > 1 and rotn:
        oscillation = (float(scan_axis[0]), float(scan_axis[1] - scan_axis[0]))
    else:
        # If not a rotation, or only one image, => stills, oscillation range = 0
        angle = float(scan_axis[0]) if rotn else 0
        oscillation = (angle, 0)

    # Get the exposure time
    if "frame_time" in detector_obj.handle:
        frame_time = float(detector_obj.handle["frame_time"][()])
        exposure_time = flex.double(num_images, frame_time)
        epochs = flex.double(num_images)
        for i in range(1, len(epochs)):
            epochs[i] = epochs[i - 1] + exposure_time[i - 1]
    else:
        exposure_time = flex.double(num_images, 0)
        epochs = flex.double(num_images, 0)

    # Construct the model
    return Scan(image_range, oscillation, exposure_time, epochs)


class CrystalFactory:
    """
    A class to create a crystal model from NXmx stuff
    """

    def __init__(self, obj):
        # Get the crystal parameters
        unit_cell_parameters = list(obj.handle["unit_cell"][0])
        unit_cell = cctbx.uctbx.unit_cell(unit_cell_parameters)
        U = list(obj.handle["orientation_matrix"][0].flatten())
        U = matrix.sqr(U)
        B = matrix.sqr(unit_cell.fractionalization_matrix()).transpose()
        A = U * B
        Ai = A.inverse()
        real_space_a = Ai[0:3]
        real_space_b = Ai[3:6]
        real_space_c = Ai[6:9]

        # Get the space group symbol
        space_group_symbol = obj.handle["unit_cell_group"][()]

        # Create the model
        self.model = Crystal(
            real_space_a, real_space_b, real_space_c, space_group_symbol
        )


class DetectorGroupDataList:
    """
    A class to make it easier to access the data from multiple datasets.
    This version brings in all the panels from a detector group with several detectors.
    """

    def __init__(self, datalists):
        self._datalists = datalists
        lengths = [len(datalist) for datalist in datalists]
        self._num_images = lengths[0]
        assert all(length == self._num_images for length in lengths), (
            "Not all datasets are the same length"
        )

    def __len__(self):
        return self._num_images

    def __getitem__(self, index):
        return tuple(itertools.chain.from_iterable(dl[index] for dl in self._datalists))


def get_detector_module_slices(detector):
    """
    Helper function to read data_origin and data_size from the NXdetector_modules in a
    NXdetector.  Returns a list of lists, where each sublist is a list of slices in
    slow to fast order.
    Assumes slices are stored in NeXus in slow to fast order.
    """
    # get the set of leaf modules (handles nested modules)
    modules = []
    for nx_detector_module in detector.modules:
        if not find_class(nx_detector_module.handle, "NXdetector_module"):
            modules.append(nx_detector_module)

    all_slices = []
    for module in modules:
        data_origin = module.handle["data_origin"]
        data_size = module.handle["data_size"]
        all_slices.append(
            [
                slice(int(start), int(start + step), 1)
                for start, step in zip(data_origin, data_size)
            ]
        )
    return all_slices


class MultiPanelDataList:
    """
    A class to make it easier to access the data from multiple datasets.
    Also handles multi-panel data as described in a series of NXdetector_modules
    """

    def __init__(self, datasets, detector):
        self._datasets = datasets
        self._num_images = 0
        self._lookup = []
        self._offset = [0]
        for i, dataset in enumerate(self._datasets):
            self._num_images += dataset.shape[0]
            self._lookup.extend([i] * dataset.shape[0])
            self._offset.append(self._num_images)

        self._all_slices = get_detector_module_slices(detector)

    def __len__(self):
        return self._num_images

    def __getitem__(self, index):
        d = self._lookup[index]
        i = index - self._offset[d]

        all_data = []

        for module_slices in self._all_slices:
            slices = [slice(i, i + 1, 1)]
            slices.extend(module_slices)
            data_as_flex = dataset_as_flex(self._datasets[d], tuple(slices))
            data_as_flex.reshape(
                flex.grid(data_as_flex.all()[-2:])
            )  # handle 3 or 4 dimension arrays
            all_data.append(data_as_flex)
        return tuple(all_data)


DataFactoryCache = collections.namedtuple(
    "DataFactoryCache", "ndim shape filename is_virtual"
)


class DataFactory:
    """
    A class to make it easier to access data from multiple datasets.
    """

    def __init__(self, obj, max_size=0, cached_information=None):
        """
        cached_information is a dictionary of DataFactoryCache named tuples.
        The dictionary key corresponds to the object handle key.
        """

        self.clear_cache()

        DataSetInformation = collections.namedtuple(
            "DataSetInformation", "accessor file shape"
        )
        datasets = []
        for key in sorted(obj.handle):
            if key.startswith("_filename_"):
                continue

            if cached_information and key in cached_information:
                ohk = cached_information[key]
                filename = ohk.filename
            else:
                ohk = obj.handle[key]  # this opens the file
                filename = ohk.file.filename

            # datasets in this context mean ones which contain diffraction images
            # so must have ndim > 1 - for example omega can also be nexus data set
            # but with 1 dimension...
            if ohk.ndim == 1:
                continue

            dsi = DataSetInformation(
                accessor=(lambda obj=obj, key=key: obj.handle[key]),
                file=filename,
                shape=ohk.shape,
            )
            if ohk.is_virtual:
                datasets = [dsi]
                break
            else:
                datasets.append(dsi)

        self._datasets = tuple(datasets)
        self._num_images = 0
        self._lookup = []
        self._offset = [0]
        self._shape = None

        if len(self._datasets) == 1 and max_size:
            self._shape = self._datasets[0].shape
            self._num_images = max_size
            self._lookup.extend([0] * max_size)
            self._offset.append(max_size)
        else:
            for i, dataset in enumerate(self._datasets):
                if self._shape is None:
                    self._shape = dataset.shape
                else:
                    assert self._shape[-2:] == dataset.shape[-2:]
                    self._shape = (self._shape[0] + dataset.shape[0],) + self._shape[1:]
                self._num_images += dataset.shape[0]
                self._lookup.extend([i] * dataset.shape[0])
                self._offset.append(self._num_images)

    def clear_cache(self):
        self._cache = (None, None)

    def __len__(self):
        return self._num_images

    def shape(self):
        return self._shape

    def __getitem__(self, index):
        d = self._lookup[index]
        i = index - self._offset[d]

        # a lock-free most-recently-used file handle cache based on
        # immutability of python tuples
        cached_handle = self._cache
        if self._datasets[d].file == cached_handle[0]:
            data = cached_handle[1]
        else:
            data = self._datasets[d].accessor()
            self._cache = (self._datasets[d].file, data)

        N, height, width = self._datasets[d].shape
        data_as_flex = dataset_as_flex(
            data, (slice(i, i + 1, 1), slice(0, height, 1), slice(0, width, 1))
        )
        data_as_flex.reshape(flex.grid(data_as_flex.all()[1:]))
        return data_as_flex


def detectorgroupdatafactory(obj, instrument):
    """Function to handle reading data from a detector with a NXdetector_group"""

    mapping = {}
    for key in sorted(obj.handle):
        if key.startswith("_filename_"):
            continue

        # When multiple datasets are in NXData, it is required to specify one as the signal
        # In such a case, we skip all datasets except the one flagged as signal
        if "signal" in obj.handle.attrs and key != obj.handle.attrs["signal"]:
            continue

        # datasets in this context mean ones which contain diffraction images
        # so must have ndim > 1 - for example omega can also be nexus data set
        # but with 1 dimension...

        dataset = obj.handle[key]  # this opens the file
        if dataset.ndim == 1:
            continue

        # Map NXdetector names to list of datasets
        dataset_name = key
        found_it = False
        for detector in instrument.detectors:
            if dataset_name in detector.handle:
                found_it = True
                detector_name = os.path.basename(detector.handle.name)
                if detector_name in mapping:
                    assert (
                        dataset_name not in mapping[detector_name]["dataset_names"]
                    ), "Dataset %s found in > 1 NXdetectors" % dataset_name
                    mapping[detector_name]["dataset_names"].add(dataset_name)
                    mapping[detector_name]["datasets"].append(dataset)
                else:
                    mapping[detector_name] = {
                        "dataset_names": {dataset_name},
                        "datasets": [dataset],
                        "detector": detector,
                    }
        assert found_it, "Couldn't match dataset %s to a NXdetector" % dataset_name

    # Create a list of multipanel datalist objects
    return DetectorGroupDataList(
        [
            MultiPanelDataList(
                detector_mapping["datasets"], detector_mapping["detector"]
            )
            for detector_mapping in mapping.values()
        ]
    )


class MaskFactory:
    """
    A class to create an object to hold the pixel mask data
    """

    def __init__(self, objects, index=None):
        def make_mask(dset, index):
            i = 0 if index is None else index
            mask = []
            for module_slices in all_slices:
                assert len(dset.shape) in [len(module_slices), len(module_slices) + 1]
                if len(dset.shape) == len(module_slices):
                    slices = []  # single image mask
                else:
                    slices = [slice(i, i + 1, 1)]  # multi-image mask
                slices.extend(module_slices)
                data_as_flex = dataset_as_flex_int(dset.id.id, tuple(slices))
                data_as_flex.reshape(
                    flex.grid(data_as_flex.all()[-2:])
                )  # handle 3 or 4 dimension arrays
                mask.append(data_as_flex == 0)
            return tuple(mask)

        self.mask = None
        for obj in objects:
            handle = obj.handle
            if "pixel_mask_applied" in handle and handle["pixel_mask_applied"]:
                if self.mask is None:
                    self.mask = []
                if "pixel_mask" in handle:
                    shape = handle["pixel_mask"].shape
                    all_slices = get_detector_module_slices(obj)
                    if len(all_slices) == 1:
                        all_slices = [[slice(0, shape[0], 1), slice(0, shape[1], 1)]]
                    self.mask.extend(list(make_mask(handle["pixel_mask"], index)))
                elif "detectorSpecific" in handle:
                    if "pixel_mask" in handle["detectorSpecific"]:
                        shape = handle["detectorSpecific"]["pixel_mask"].shape
                        all_slices = get_detector_module_slices(obj)
                        if len(all_slices) == 1:
                            all_slices = [
                                [slice(0, shape[0], 1), slice(0, shape[1], 1)]
                            ]
                        self.mask.extend(
                            list(
                                make_mask(
                                    handle["detectorSpecific"]["pixel_mask"], index
                                )
                            )
                        )
        if self.mask is not None:
            self.mask = tuple(self.mask)
