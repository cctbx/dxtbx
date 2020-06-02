from __future__ import absolute_import, division, print_function

import collections
import itertools
import math
import os
from builtins import range

import h5py
import numpy
import six

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
)

try:
    from dxtbx_format_nexus_ext import (
        dataset_as_flex_int,
        dataset_as_flex_double,
        dataset_as_flex_float,
    )
except ImportError:
    # Workaround for psana build, which doesn't link HDF5 properly
    if "SIT_ROOT" not in os.environ:
        raise


def dataset_as_flex(dataset, selection):
    if numpy.issubdtype(dataset.dtype, numpy.integer):
        return dataset_as_flex_int(dataset.id.id, selection)
    else:
        assert numpy.issubdtype(dataset.dtype, numpy.floating)
        if dataset.dtype in [
            numpy.half,
            numpy.single,
            numpy.float_,
            numpy.float16,
            numpy.float32,
        ]:
            return dataset_as_flex_float(dataset.id.id, selection)
        elif dataset.dtype in [
            numpy.double,
            numpy.longfloat,
            numpy.float64,
            numpy.float96,
            numpy.float128,
        ]:
            return dataset_as_flex_double(dataset.id.id, selection)
        else:
            assert False, "unknown floating data type (%s)" % str(dataset.dtype)


class NXValidationError(RuntimeError):
    """A specific exception to record validation errors"""


def local_visit(nx_file, visitor):
    """Implementation of visitor to replace node.visititems

    Will dereference soft links but should avoid walking down into external files,
    I think - not a property that NXgroup.visititems(visitor) has.
    https://github.com/cctbx/dxtbx/issues/74

    Args:
      nx_file: hdf5 file node
      visitor: visitor function to act on children
    """
    for key in nx_file.keys():
        if isinstance(nx_file.get(key, getlink=True), h5py.ExternalLink):
            # Follow links but do not recurse into external files
            continue

        # Do not iterate over .values().
        # As .values() is not a true generator that would mean that all value objects
        # would have to be represented in memory at the same time. If the value
        # objects refer to external files then those are kept open until the loop
        # terminates, at which point all of the file handles are garbage collected
        # and closed at once.
        k = nx_file[key]

        if "NX_class" not in k.attrs:
            continue
        visitor(k.name, k)
        local_visit(k, visitor)


def find_entries(nx_file, entry):
    """
    Find NXmx entries
    """
    hits = []

    def visitor(name, obj):
        if "NX_class" in obj.attrs:
            if numpy.string_(obj.attrs["NX_class"]) in [
                numpy.string_("NXentry"),
                numpy.string_("NXsubentry"),
            ]:
                if "definition" in obj:
                    if obj["definition"][()] == numpy.string_("NXmx"):
                        hits.append(obj)

    visitor(entry, nx_file[entry])
    local_visit(nx_file, visitor)
    return hits


def find_class(nx_file, nx_class):
    """
    Find a given NXclass
    """
    hits = []
    nx_class = numpy.string_(nx_class)

    def visitor(name, obj):
        if numpy.string_("NX_class") in obj.attrs:
            if numpy.string_(obj.attrs["NX_class"]) == nx_class:
                hits.append(obj)

    local_visit(nx_file, visitor)
    return hits


def convert_units(value, input_units, output_units):
    """
    Hacky utility function to convert units
    """
    if six.PY3 and isinstance(input_units, bytes):
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
        raise RuntimeError("Can't convert units %r to %r" % (input_units, output_units))


def visit_dependencies(nx_file, item, visitor=None):
    """
    Walk the dependency chain and call a visitor function
    """
    dependency_chain = set()
    if os.path.basename(item) == "depends_on":
        depends_on = nx_file[item][()]
    else:
        depends_on = nx_file[item].attrs["depends_on"]
    while not depends_on == numpy.string_("."):
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
            depends_on = nx_file[depends_on].attrs["depends_on"]
        except Exception:
            raise RuntimeError("'%s' contains no depends_on attribute" % depends_on)


def construct_vector(nx_file, item, vector=None):
    """
    Walk the dependency chain and create the absolute vector
    """

    class TransformVisitor(object):
        def __init__(self, vector):
            self.vector = matrix.col(vector)

        def visit(self, nx_file, depends_on):
            item = nx_file[depends_on]
            value = item[()]
            units = item.attrs["units"]
            ttype = item.attrs["transformation_type"]
            vector = matrix.col(item.attrs["vector"])
            if ttype == numpy.string_("translation"):
                value = convert_units(value, units, "mm")
                if hasattr(value, "__iter__") and len(value) == 1:
                    value = value[0]
                self.vector = vector * value + self.vector
            elif ttype == numpy.string_("rotation"):
                if hasattr(value, "__iter__") and len(value):
                    value = value[0]
                if units == "rad":
                    deg = False
                elif units == "deg":
                    deg = True
                else:
                    raise RuntimeError("Invalid units: %s" % units)
                self.vector.rotate(axis=vector, angle=value, deg=deg)
            else:
                raise RuntimeError("Unknown transformation_type: %s" % ttype)

    if vector is None:
        value = nx_file[item][()]
        units = nx_file[item].attrs["units"]
        ttype = nx_file[item].attrs["transformation_type"]
        vector = nx_file[item].attrs["vector"]
        if "offset" in nx_file[item].attrs:
            offset = nx_file[item].attrs["offset"]
            offset = convert_units(offset, units, "mm")
        else:
            offset = vector * 0.0
        if ttype == numpy.string_("translation"):
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

    class Visitor(object):
        def __init__(self):
            self._axes = flex.vec3_double()
            self._angles = flex.double()
            self._axis_names = flex.std_string()
            self._is_scan_axis = flex.bool()

        def visit(self, nx_file, depends_on):
            item = nx_file[depends_on]
            value = item[()]
            units = item.attrs["units"]
            ttype = item.attrs["transformation_type"]
            vector = [float(v) for v in item.attrs["vector"]]
            if ttype == numpy.string_("translation"):
                return
            elif ttype == numpy.string_("rotation"):
                if hasattr(value, "__iter__") and len(value):
                    value = value[0]
                if units == numpy.string_("rad"):
                    value *= math.pi / 180
                elif units not in [
                    numpy.string_("deg"),
                    numpy.string_("degree"),
                    numpy.string_("degrees"),
                ]:
                    raise RuntimeError("Invalid units: %s" % units)

                # is the axis moving? Check the values for this axis
                v = item[()]
                if min(v) < max(v):
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
                assert (
                    self._is_scan_axis.count(True) == 1
                ), "Only one axis can be a scan axis: %s" % list(self._is_scan_axis)
                scan_axis = flex.first_index(self._is_scan_axis, True)

            # Rotate 180 about up from McStas coordinate system
            cb_op = (-1, 0, 0, 0, 1, 0, 0, 0, -1)
            return cb_op * self._axes, self._angles, self._axis_names, scan_axis

    if vector is None:
        value = nx_file[item][()]
        units = nx_file[item].attrs["units"]
        ttype = nx_file[item].attrs["transformation_type"]
        vector = nx_file[item].attrs["vector"]
        if "offset" in nx_file[item].attrs:
            offset = nx_file[item].attrs["offset"]
        else:
            offset = vector * 0.0
        if ttype == numpy.string_("translation"):
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


class NXdetector_module(object):
    """
    A class to hold a handle to NXdetector_module
    """

    def __init__(self, handle):
        self.handle = handle


class NXdetector_group(object):
    """
    A class to hold a handle to NXdetector_group
    """

    def __init__(self, handle):
        self.handle = handle


class NXdetector(object):
    """
    A class to handle a handle to NXdetector
    """

    def __init__(self, handle):
        self.handle = handle

        # Find the NXdetector_modules
        self.modules = []
        for entry in find_class(self.handle, "NXdetector_module"):
            self.modules.append(NXdetector_module(entry))

        # Check we've got some stuff
        if not self.modules:
            raise NXValidationError("No NXdetector_module in %s" % self.handle.name)


class NXinstrument(object):
    """
    A class to hold a handle to NXinstrument
    """

    def __init__(self, handle):
        self.handle = handle

        # Find the NXdetector
        self.detectors = []
        for entry in find_class(self.handle, "NXdetector"):
            self.detectors.append(NXdetector(entry))

        # Check we've got stuff
        if not self.detectors:
            raise NXValidationError("No NXdetector in %s" % self.handle.name)

        # Find any detector groups
        self.detector_groups = []
        for entry in find_class(self.handle, "NXdetector_group"):
            self.detector_groups.append(NXdetector_group(entry))

        # Find the NXbeam
        self.beams = []
        for entry in find_class(self.handle, "NXbeam"):
            self.beams.append(NXbeam(entry))


class NXbeam(object):
    """
    A class to hold a handle to NXbeam
    """

    def __init__(self, handle):
        self.handle = handle


class NXsample(object):
    """
    A class to hold a handle to NXsample
    """

    def __init__(self, handle):
        self.handle = handle

        # Find the NXbeam
        self.beams = []
        for entry in find_class(self.handle, "NXbeam"):
            self.beams.append(NXbeam(entry))


class NXdata(object):
    """
    A class to hold a handle to NXdata
    """

    def __init__(self, handle):
        self.handle = handle


class NXmxEntry(object):
    """
    A class to hold a handle to NXmx entries
    """

    def __init__(self, handle):
        self.handle = handle

        # Find the NXinstrument
        self.instruments = []
        for entry in find_class(self.handle, "NXinstrument"):
            self.instruments.append(NXinstrument(entry))

        # Find the NXsample
        self.samples = []
        for entry in find_class(self.handle, "NXsample"):
            self.samples.append(NXsample(entry))

        # Find the NXdata
        self.data = []
        for entry in find_class(self.handle, "NXdata"):
            self.data.append(NXdata(entry))

        # Check we've got some stuff
        if not self.instruments:
            raise NXValidationError("No NXinstrument in %s" % self.handle.name)
        if not self.samples:
            raise NXValidationError("No NXsample in %s" % self.handle.name)
        if not self.data:
            raise NXValidationError("No NXdata in %s" % self.handle.name)


class NXmxReader(object):
    """
    A hacky class to read an NXmx file
    """

    def __init__(self, filename=None, handle=None):
        # Get the file handle
        if filename is not None:
            handle = h5py.File(filename, "r")

        # Find the NXmx entries
        self.entries = []
        for entry in find_entries(handle, "/"):
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
        return bool(find_entries(handle, "/"))


class BeamFactory(object):
    """
    A class to create a beam model from NXmx stuff
    """

    def __init__(self, obj, index=None):
        # Get the items from the NXbeam class
        wavelength = obj.handle["incident_wavelength"]
        wavelength_value = wavelength[()]
        wavelength_units = wavelength.attrs["units"]

        if (
            index is not None
            and hasattr(wavelength_value, "__iter__")
            and len(wavelength_value) > 1
        ):
            wavelength_value = wavelength_value[index]

        # Convert wavelength to Angstroms
        wavelength_value = float(
            convert_units(wavelength_value, wavelength_units, "angstrom")
        )

        # Construct the beam model
        self.model = Beam(direction=(0, 0, 1), wavelength=wavelength_value)


class BeamFactory_(object):
    """
    A class to create a beam model from NXmx stuff, backported from DIALS 3.0
    """

    def __init__(self, obj):
        self.obj = obj
        self.model = None
        self.index = None

    def load_model(self, index=None):
        # Cached model
        if self.model is not None and index == self.index:
            return self.model

        # Get the items from the NXbeam class
        wavelength = self.obj.handle["incident_wavelength"]
        wavelength_weights = self.obj.handle.get("incident_wavelength_weights")
        if wavelength.shape in (tuple(), (1,)):
            wavelength_value = wavelength[()]
        elif len(wavelength.shape) == 1:
            if wavelength_weights is None:
                if index is None:
                    index = 0
                wavelength_value = wavelength[index]
            else:
                raise NotImplementedError("Spectra not implemented")
        else:
            raise NotImplementedError("Spectra not implemented")
        wavelength_units = wavelength.attrs["units"]

        # Convert wavelength to Angstroms
        wavelength_value = float(
            convert_units(wavelength_value, wavelength_units, "angstrom")
        )

        # Construct the beam model
        self.index = index
        self.model = Beam(direction=(0, 0, 1), wavelength=wavelength_value)
        return self.model


def get_change_of_basis(transformation):
    """
    Get the 4x4 homogenous coordinate matrix for a given NXtransformation.
    """
    # Change of basis to convert from NeXus to IUCr/ImageCIF convention
    n2i_cob = sqr((-1, 0, 0, 0, 1, 0, 0, 0, -1))

    axis_type = numpy.string_(transformation.attrs["transformation_type"])

    vector = n2i_cob * col(transformation.attrs["vector"]).normalize()
    setting = transformation[0]
    units = numpy.string_(transformation.attrs["units"])

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

    if axis_type == numpy.string_("rotation"):
        if units == numpy.string_("rad"):
            deg = False
        elif units in [
            numpy.string_("deg"),
            numpy.string_("degree"),
            numpy.string_("degrees"),
        ]:
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
    elif axis_type == numpy.string_("translation"):
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
        raise ValueError("Unrecognized tranformation type: %s" % axis_type)

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
        parent_id = numpy.string_(current.attrs["depends_on"])

        if parent_id == numpy.string_("."):
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
    transformation's depends_on is ".".  Parent is the tranformation that the
    top level transformation in this chain of transformations depends on.
    """

    cob = get_change_of_basis(transformation)

    parent_id = numpy.string_(transformation.attrs["depends_on"])

    if parent_id == numpy.string_("."):
        return None, cob
    parent = transformation.parent[parent_id]

    if "equipment_component" in transformation.attrs:
        eq_comp = transformation.attrs["equipment_component"]
        parent_eq_comp = parent.attrs["equipment_component"]
        if eq_comp == parent_eq_comp:
            non_matching_parent, parent_cob = get_cumulative_change_of_basis(parent)
            return non_matching_parent, parent_cob * cob

    return parent, cob


class DetectorFactoryFromGroup(object):
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
            ], "Hierarchy of detectors not supported. Hierarchy of module components within detector elements is supported"

            if parent_id == -1:
                assert root_name is None, "Multiple roots not supported"
                root_name = group_names[i]
            else:
                expected_detectors.append(group_names[i].astype(str))

        assert root_name is not None, "Detector root not found"
        assert sorted(
            os.path.basename(d.handle.name) for d in instrument.detectors
        ) == sorted(
            expected_detectors
        ), "Mismatch between detector group names and detectors available"

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
            modules = []
            for nx_detector_module in nx_detector.modules:
                if not find_class(nx_detector_module.handle, "NXdetector_module"):
                    modules.append(nx_detector_module)

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

                # Get the trusted range of pixel values
                underload = (
                    float(nx_detector.handle["underload_value"][()])
                    if "underload_value" in nx_detector.handle
                    else -400
                )
                overload = (
                    float(nx_detector.handle["saturation_value"][()])
                    if "saturation_value" in nx_detector.handle
                    else 90000
                )
                trusted_range = underload, overload

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
                pg = None
                for transform in reversed(chain):
                    name = str(os.path.basename(transform.name))
                    if (
                        pg is None
                    ):  # The first transform will be the root of the hiearchy
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
                parent, cob = get_cumulative_change_of_basis(depends_on)
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
                    material = str(nx_detector.handle["sensor_material"][()])
                    p.set_material(material)

                    # Compute the attenuation coefficient.
                    # This will fail for undefined composite materials
                    # mu_at_angstrom returns cm^-1, but need mu in mm^-1
                    if material == "Si":
                        pass
                    elif material == "Silicon":
                        material = "Si"
                    elif material == "Sillicon":
                        material = "Si"
                    elif material == "CdTe":
                        pass
                    elif material == "GaAs":
                        pass
                    else:
                        raise RuntimeError("Unknown material: %s" % material)
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


def known_backwards(image_size):
    """
    Tests for special cases for known data where image size is backwards from NeXus spec

    image_size is in dxtbx image size order (fast, slow)
    """
    return image_size in [
        (4362, 4148),  # Eiger 2X 16M @ DLS
        (4371, 4150),  # Eiger 16M @ Spring8
        (2162, 2068),  # Eiger 2x 4M @ VMXi, DLS
        (2167, 2070),  # Eiger 1 4M @ VMXi, DLS
        (3269, 3110),  # Eiger 9M Proxima2A beamline @ SOLEIL
    ]


class DetectorFactory(object):
    """
    A class to create a detector model from NXmx stuff
    """

    def __init__(self, obj, beam):
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

        # Get the trusted range of pixel values
        if "saturation_value" in nx_detector:
            trusted_range = (-1, float(nx_detector["saturation_value"][()]))
        else:
            trusted_range = (-1, 99999999)

        # Get the detector thickness
        thickness = nx_detector["sensor_thickness"]
        thickness_value = float(thickness[()])
        thickness_units = thickness.attrs["units"]
        thickness_value = float(convert_units(thickness_value, thickness_units, "mm"))

        # Get the detector material
        material = {
            numpy.string_("Si"): "Si",
            numpy.string_("Silicon"): "Si",
            numpy.string_("Sillicon"): "Si",
            numpy.string_("CdTe"): "CdTe",
            numpy.string_("GaAs"): "GaAs",
        }.get(nx_detector["sensor_material"][()])
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
        mu = table.mu_at_angstrom(wavelength) / 10.0

        # Construct the detector model
        pixel_size = (fast_pixel_direction_value, slow_pixel_direction_value)
        # image size stored slow to fast but dxtbx needs fast to slow
        image_size = tuple(int(x) for x in reversed(nx_module["data_size"][-2:]))

        if known_backwards(image_size):
            image_size = tuple(reversed(image_size))

        self.model = Detector()
        self.model.add_panel(
            Panel(
                detector_type,
                detector_name,
                tuple(fast_axis),
                tuple(slow_axis),
                tuple(origin),
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


class GoniometerFactory(object):
    """
    A class to create a goniometer model from NXmx stuff
    """

    def __init__(self, obj):
        if obj.handle["depends_on"][()] == ".":
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
    if obj.handle["depends_on"][()] == ".":
        return
    thing = obj.handle.file[obj.handle["depends_on"][()]]
    tree = get_depends_on_chain_using_equipment_components(thing)
    for t in tree:
        o = obj.handle.file[t.name]
        if o.attrs["transformation_type"] == numpy.string_("rotation"):
            # if this is changing, assume is scan axis
            v = o[()]
            if min(v) < max(v):
                return o
    raise ValueError("no rotation found")


def find_scanning_axis(obj):
    if obj.handle["depends_on"][()] == ".":
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
    if obj.handle["depends_on"][()] == ".":
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

    rotn = scan_axis.attrs["transformation_type"] == numpy.string_("rotation")

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


class ScanFactory(object):
    """
    A class to create a scan model from NXmx stuff
    """

    def __init__(self, obj, detector_obj):
        # Get the image and oscillation range - need to search for rotations
        # in dependency tree - if not, find translations or just the thing
        # the sample depends on
        try:
            scan_axis = find_goniometer_rotation(obj)
        except ValueError:
            scan_axis = find_scanning_axis(obj)

        if scan_axis is None:
            scan_axis = obj.handle.file[obj.handle["depends_on"][()]]

        image_range = (1, len(scan_axis))

        rotn = scan_axis.attrs["transformation_type"] == numpy.string_("rotation")

        if len(scan_axis) > 1:
            if rotn:
                oscillation = (float(scan_axis[0]), float(scan_axis[1] - scan_axis[0]))
            else:
                # nothing to see here - should really work out how to dig out
                # the setting for this - but not today
                oscillation = (0, 0)
            is_sequence = True
        else:
            oscillation = (float(scan_axis[0]), 0.0)
            is_sequence = False

        # Get the exposure time
        num_images = len(scan_axis)
        if "frame_time" in detector_obj.handle:
            frame_time = float(detector_obj.handle["frame_time"][()])
            exposure_time = flex.double(num_images, frame_time)
            epochs = flex.double(num_images)
            for i in range(1, len(epochs)):
                epochs[i] = epochs[i - 1] + exposure_time[i - 1]
        else:
            exposure_time = flex.double(num_images, 0)
            epochs = flex.double(num_images, 0)

        if is_sequence is True:

            # Construct the model
            self.model = Scan(image_range, oscillation, exposure_time, epochs)

        else:

            # if this is not a sequence, then this is stills, in which case... should
            # we have scans? and, even if we do, we have a small problem in that
            # this returns a list of scans which even less works...

            self.model = []
            for i, image in enumerate(range(image_range[0], image_range[1] + 1)):
                self.model.append(
                    Scan(
                        (image, image),
                        oscillation,
                        exposure_time[i : i + 1],
                        epochs[i : i + 1],
                    )
                )

            self.model = None


class CrystalFactory(object):
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


class DetectorGroupDataList(object):
    """
    A class to make it easier to access the data from multiple datasets.
    This version brings in all the panels from a detector group with several detectors.
    """

    def __init__(self, datalists):
        self._datalists = datalists
        lengths = [len(datalist) for datalist in datalists]
        self._num_images = lengths[0]
        assert all(
            length == self._num_images for length in lengths
        ), "Not all datasets are the same length"

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
        if known_backwards((data_size[-1], data_size[-2])):
            all_slices.append(
                list(
                    reversed(
                        [
                            slice(int(start), int(start + step), 1)
                            for start, step in zip(data_origin, data_size)
                        ]
                    )
                )
            )
        else:
            all_slices.append(
                [
                    slice(int(start), int(start + step), 1)
                    for start, step in zip(data_origin, data_size)
                ]
            )
    return all_slices


class MultiPanelDataList(object):
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


DataFactoryCache = collections.namedtuple("DataFactoryCache", "ndim shape filename")


class DataFactory(object):
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

            datasets.append(
                DataSetInformation(
                    accessor=(lambda obj=obj, key=key: obj.handle[key]),
                    file=filename,
                    shape=ohk.shape,
                )
            )

        self._datasets = tuple(datasets)
        self._num_images = 0
        self._lookup = []
        self._offset = [0]

        if len(self._datasets) == 1 and max_size:
            self._num_images = max_size
            self._lookup.extend([0] * max_size)
            self._offset.append(max_size)
        else:
            for i, dataset in enumerate(self._datasets):
                self._num_images += dataset.shape[0]
                self._lookup.extend([i] * dataset.shape[0])
                self._offset.append(self._num_images)

    def clear_cache(self):
        self._cache = (None, None)

    def __len__(self):
        return self._num_images

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

        # datasets in this context mean ones which contain diffraction images
        # so must have ndim > 1 - for example omega can also be nexus data set
        # but with 1 dimension...

        dataset = obj.handle[key]  # this opens the file
        if dataset.ndim == 1:
            continue

        # Map NXdetector names to list of datasets
        dataset_name = os.path.basename(dataset.name)
        found_it = False
        for detector in instrument.detectors:
            if dataset_name in detector.handle:
                found_it = True
                detector_name = os.path.basename(detector.handle.name)
                if detector_name in mapping:
                    assert (
                        dataset_name not in mapping[detector_name]["dataset_names"]
                    ), ("Dataset %s found in > 1 NXdetectors" % dataset_name)
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


class MaskFactory(object):
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
                    all_slices = get_detector_module_slices(obj)
                    self.mask.extend(list(make_mask(handle["pixel_mask"], index)))
                elif "detectorSpecific" in handle:
                    if "pixel_mask" in handle["detectorSpecific"]:
                        all_slices = get_detector_module_slices(obj)
                        self.mask.extend(
                            list(
                                make_mask(
                                    handle["detectorSpecific"]["pixel_mask"], index
                                )
                            )
                        )
        if self.mask is not None:
            self.mask = tuple(self.mask)
