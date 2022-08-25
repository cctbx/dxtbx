from __future__ import annotations

import dataclasses
import datetime
import logging
import operator
from collections import abc, namedtuple
from functools import cached_property, reduce
from typing import Iterable, Iterator, Sequence, Union, overload

import dateutil.parser
import h5py
import numpy as np
import pint
from scipy.spatial.transform import Rotation

# NeXus field type for type annotations
# https://manual.nexusformat.org/nxdl-types.html#nxdl-field-types-and-units
NXBoolT = Union[bool, np.ndarray]
NXFloatT = Union[float, np.ndarray]
NXIntT = Union[int, np.ndarray]
NXNumberT = Union[NXFloatT, NXIntT]


ureg = pint.UnitRegistry()


logger = logging.getLogger(__name__)


NXNode = Union[h5py.File, h5py.Group]


class NXNumber(abc.Sequence):
    def __init__(self, handle: h5py.Dataset, unit: pint.Unit | None):
        self._handle = handle
        self._unit = unit

    def __getitem__(self, key) -> NXNumberT:
        if self._unit:
            return self._handle[key] * self._unit
        return self._handle[key]

    def __len__(self):
        return len(self._handle)


def h5str(h5_value: str | np.bytes_ | bytes | None) -> str | None:
    """
    Convert a value returned from an h5py attribute to str.

    h5py can return either a bytes-like (numpy.string_) or str object
    for attribute values depending on whether the value was written as
    fixed or variable length. This function collapses the two to str.
    """
    if isinstance(h5_value, (np.bytes_, bytes)):
        return h5_value.decode("utf-8")
    return h5_value


def units(data: h5py.Dataset, default: str | None = None) -> pint.Unit:
    """Extract the units attribute, if any, from an h5py data set."""
    return ureg.Unit(h5str(data.attrs.get("units", default)))


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


class H5Mapping(abc.Mapping):
    def __init__(self, handle: h5py.File | h5py.Group):
        self._handle = handle

    def __getitem__(self, key: str) -> h5py.Group | h5py.Dataset:
        return self._handle[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._handle)

    def __len__(self) -> int:
        return len(self._handle)

    @cached_property
    def path(self) -> str | None:
        return h5str(self._handle.name)


class NXmx(H5Mapping):
    def __init__(self, handle):
        super().__init__(handle)
        self._entries = [
            entry
            for entry in find_class(handle, "NXentry")
            if "definition" in entry and h5str(entry["definition"][()]) == "NXmx"
        ]

    @cached_property
    def entries(self) -> list[NXentry]:
        return [NXentry(entry) for entry in self._entries]


class NXentry(H5Mapping):
    """NXentry describes the measurement.

    The top-level NeXus group which contains all the data and associated information
    that comprise a single measurement. It is mandatory that there is at least one group
    of this type in the NeXus file.
    """

    def __init__(self, handle):
        super().__init__(handle)
        self._data, self._instruments, self._samples, self._sources = find_classes(
            handle, "NXdata", "NXinstrument", "NXsample", "NXsource"
        )

    @cached_property
    def instruments(self) -> list[NXinstrument]:
        return [NXinstrument(instrument) for instrument in self._instruments]

    @cached_property
    def samples(self) -> list[NXsample]:
        return [NXsample(sample) for sample in self._samples]

    @cached_property
    def data(self) -> list[NXdata]:
        return [NXdata(data) for data in self._data]

    @cached_property
    def source(self) -> NXsource:
        return NXsource(self._sources[0])

    @cached_property
    def start_time(self) -> datetime.datetime:
        """Starting time of measurement.

        ISO 8601 time/date of the first data point collected in UTC, using the Z suffix
        to avoid confusion with local time. Note that the time zone of the beamline
        should be provided in NXentry/NXinstrument/time_zone.
        """
        if "start_time" in self._handle:
            return dateutil.parser.isoparse(h5str(self._handle["start_time"][()]))

    @cached_property
    def end_time(self) -> datetime.datetime | None:
        """Ending time of measurement.

        ISO 8601 time/date of the last data point collected in UTC, using the Z suffix
        to avoid confusion with local time. Note that the time zone of the beamline
        should be provided in NXentry/NXinstrument/time_zone. This field should only be
        filled when the value is accurately observed. If the data collection aborts or
        otherwise prevents accurate recording of the end_time, this field should be
        omitted.
        """
        if "end_time" in self._handle:
            return dateutil.parser.isoparse(h5str(self._handle["end_time"][()]))
        return None

    @cached_property
    def end_time_estimated(self) -> datetime.datetime:
        """Estimated ending time of the measurement.

        ISO 8601 time/date of the last data point collected in UTC, using the Z suffix
        to avoid confusion with local time. Note that the time zone of the beamline
        should be provided in NXentry/NXinstrument/time_zone. This field may be filled
        with a value estimated before an observed value is available.
        """
        if "end_time_estimated" in self._handle:
            return dateutil.parser.isoparse(
                h5str(self._handle["end_time_estimated"][()])
            )

    @cached_property
    def definition(self) -> str:
        """NeXus NXDL schema to which this file conforms."""
        return h5str(self._handle["definition"][()])


class NXdata(H5Mapping):
    """NXdata describes the plottable data and related dimension scales."""

    @cached_property
    def signal(self) -> str | None:
        """Declares which dataset is the default.

        The value is the name of the dataset to be plotted. A field of this name must
        exist (either as dataset or as a link to a dataset).

        It is recommended (as of NIAC2014) to use this attribute rather than adding a
        signal attribute to the dataset. See
        https://www.nexusformat.org/2014_How_to_find_default_data.html for a summary of
        the discussion.
        """
        return self._handle.attrs.get("signal")


class NXtransformations(H5Mapping):
    """Collection of axis-based translations and rotations to describe a geometry.

    May also contain axes that do not move and therefore do not have a transformation
    type specified, but are useful in understanding coordinate frames within which
    transformations are done, or in documenting important directions, such as the
    direction of gravity.

    A nested sequence of transformations lists the translation and rotation steps needed
    to describe the position and orientation of any movable or fixed device.

    There will be one or more transformations (axes) defined by one or more fields for
    each transformation. The all-caps name AXISNAME designates the particular axis
    generating a transformation (e.g. a rotation axis or a translation axis or a general
    axis). The attribute units="NX_TRANSFORMATION" designates the units will be
    appropriate to the transformation_type attribute:
      - NX_LENGTH for translation
      - NX_ANGLE for rotation
      - NX_UNITLESS for axes for which no transformation type is specified

    This class will usually contain all axes of a sample stage or goniometer or a
    detector. The NeXus default McSTAS coordinate frame is assumed, but additional
    useful coordinate axes may be defined by using axes for which no transformation type
    has been specified.

    The entry point (depends_on) will be outside of this class and point to a field in
    here. Following the chain may also require following depends_on links to
    transformations outside, for example to a common base table. If a relative path is
    given, it is relative to the group enclosing the depends_on specification.

    For a chain of three transformations, where T1 depends on T2  and that in turn
    depends on T3, the final transformation Tf is
        Tf = T3 T2 T1

    In explicit terms, the transformations are a subset of affine transformations
    expressed as 4x4 matrices that act on homogeneous coordinates, w = (x, y, z, 1)^T.

    For rotation and translation,
        Tr = (R  o)
             (03 1)
        Tr = (I3 t + o)
             (03 1)

    where R is the usual 3x3 rotation matrix, o is an offset vector, 03 is a row of 3
    zeros, I3 is the 3x3 identity matrix and t is the translation vector.

    o is given by the offset attribute, t is given by the vector attribute multiplied by
    the field value, and R is defined as a rotation about an axis in the direction of
    vector, of angle of the field value.

    NOTE:

    One possible use of NXtransformations is to define the motors and transformations
    for a diffractometer (goniometer). Such use is mentioned in the NXinstrument base
    class. Use one NXtransformations group for each diffractometer and name the group
    appropriate to the device. Collecting the motors of a sample table or xyz-stage in
    an NXtransformation group is equally possible.
    """

    def __init__(self, handle):
        super().__init__(handle)
        self._axes = {
            k: NXtransformationsAxis(v)
            for k, v in handle.items()
            if isinstance(v, h5py.Dataset) and "vector" in v.attrs
        }

    @cached_property
    def default(self) -> str:
        """Declares which child group contains a path leading to a NXdata group.

        It is recommended (as of NIAC2014) to use this attribute to help define the path
        to the default dataset to be plotted. See
        https://www.nexusformat.org/2014_How_to_find_default_data.html for a summary of
        the discussion.
        """
        return h5str(self._handle.attrs.get("default"))

    @cached_property
    def axes(self) -> dict[str, NXtransformationsAxis]:
        return self._axes


class NXtransformationsAxis:
    """Axis-based translation and rotation to describe a given transformation.

    For a chain of three transformations, where T1 depends on T2  and that in turn
    depends on T3, the final transformation Tf is
        Tf = T3 T2 T1

    In explicit terms, the transformations are a subset of affine transformations
    expressed as 4x4 matrices that act on homogeneous coordinates, w = (x, y, z, 1)^T.

    For rotation and translation,
        Tr = (R  o)
             (03 1)
        Tr = (I3 t + o)
             (03 1)

    where R is the usual 3x3 rotation matrix, o is an offset vector, 03 is a row of 3
    zeros, I3 is the 3x3 identity matrix and t is the translation vector.

    o is given by the offset attribute, t is given by the vector attribute multiplied by
    the field value, and R is defined as a rotation about an axis in the direction of
    vector, of angle of the field value.

    NOTE:

    One possible use of NXtransformations is to define the motors and transformations
    for a diffractometer (goniometer). Such use is mentioned in the NXinstrument base
    class. Use one NXtransformations group for each diffractometer and name the group
    appropriate to the device. Collecting the motors of a sample table or xyz-stage in
    an NXtransformation group is equally possible.
    """

    def __init__(self, handle):
        self._handle = handle

    def __len__(self) -> int:
        return self._handle.size

    @cached_property
    def path(self) -> str | None:
        return h5str(self._handle.name)

    @cached_property
    def units(self) -> pint.Unit:
        """Units of the specified transformation.

        Could be any of these: NX_LENGTH, NX_ANGLE, or NX_UNITLESS

        There will be one or more transformations defined by one or more fields for each
        transformation. The units type NX_TRANSFORMATION designates the particular axis
        generating a transformation (e.g. a rotation axis or a translation axis or a
        general axis). NX_TRANSFORMATION designates the units will be appropriate to the
        type of transformation, indicated in the NXtransformations base class by the
        transformation_type value:
          - NX_LENGTH for translation
          - NX_ANGLE for rotation
          - NX_UNITLESS for axes for which no transformation type is specified.
        """
        return units(self._handle)

    @cached_property
    def transformation_type(self) -> str:
        """The type of the transformation, either translation or rotation.

        The transformation_type may be translation, in which case the values are linear
        displacements along the axis, rotation, in which case the values are angular
        rotations around the axis.

        If this attribute is omitted, this is an axis for which there is no motion to be
        specifies, such as the direction of gravity, or the direction to the source, or
        a basis vector of a coordinate frame.

        Any of these values: translation | rotation.
        """
        return h5str(self._handle.attrs.get("transformation_type"))

    @cached_property
    def vector(self) -> NXNumberT:
        """Three values that define the axis for this transformation.

        The axis should be normalized to unit length, making it dimensionless. For
        rotation axes, the direction should be chosen for a right-handed rotation with
        increasing angle. For translation axes the direction should be chosen for
        increasing displacement. For general axes, an appropriate direction should be
        chosen.
        """
        return self._handle.attrs.get("vector")

    @cached_property
    def offset(self) -> pint.Quantity | None:
        """A fixed offset applied before the transformation (three vector components).

        This is not intended to be a substitute for a fixed translation axis but, for
        example, as the mechanical offset from mounting the axis to its dependency.
        """
        if "offset" in self._handle.attrs:
            return self._handle.attrs["offset"] * self.offset_units
        return None

    @cached_property
    def offset_units(self) -> pint.Unit:
        """Units of the offset. Values should be consistent with NX_LENGTH."""
        if "offset_units" in self._handle.attrs:
            return ureg.Unit(h5str(self._handle.attrs["offset_units"]))
        # This shouldn't be the case, but DLS EIGER NeXus files include offset without
        # accompanying offset_units, so use units instead (which should strictly only
        # apply to vector, not offset).
        # See also https://jira.diamond.ac.uk/browse/MXGDA-3668
        return self.units

    @cached_property
    def depends_on(self) -> NXtransformationsAxis | None:
        depends_on = h5str(self._handle.attrs.get("depends_on"))
        if depends_on and depends_on != ".":
            return NXtransformationsAxis(self._handle.parent[depends_on])
        return None

    def __getitem__(self, key) -> pint.Quantity:
        return np.atleast_1d(self._handle)[key] * self.units

    @cached_property
    def end(self) -> NXNumber | None:
        end_name = self._handle.name + "_end"
        if end_name in self._handle.parent:
            return NXNumber(self._handle.parent[end_name], self.units)
        return None

    @cached_property
    def increment_set(self) -> pint.Quantity | None:
        increment_set_name = self._handle.name + "_increment_set"
        if increment_set_name in self._handle.parent:
            return self._handle.parent[increment_set_name][()] * self.units
        return None

    @cached_property
    def matrix(self) -> np.ndarray:
        values = self[()]
        if np.any(values):
            values = (
                values.to("mm").magnitude
                if self.transformation_type == "translation"
                else values.to("rad").magnitude
            )
        else:
            values = values.magnitude

        if self.transformation_type == "rotation":
            R = Rotation.from_rotvec(values[:, np.newaxis] * self.vector).as_matrix()
            T = np.zeros((values.size, 3))
        else:
            R = np.identity(3)
            T = values[:, np.newaxis] * self.vector

        if self.offset is not None and np.any(self.offset):
            T += self.offset.to("mm").magnitude

        A = np.repeat(np.identity(4).reshape((1, 4, 4)), values.size, axis=0)
        A[:, :3, :3] = R
        A[:, :3, 3] = T
        return A


class NXsample(H5Mapping):
    """Any information on the sample.

    This could include scanned variables that are associated with one of the data
    dimensions, e.g. the magnetic field, or logged data, e.g. monitored temperature vs
    elapsed time.
    """

    def __init__(self, handle):
        super().__init__(handle)
        self._transformations = find_class(handle, "NXtransformations")

    @cached_property
    def name(self) -> str:
        """Descriptive name of sample"""
        return h5str(self._handle["name"][()])

    @cached_property
    def depends_on(self) -> NXtransformationsAxis | None:
        """The axis on which the sample position depends"""
        if "depends_on" in self._handle:
            depends_on = h5str(self._handle["depends_on"][()])
            if depends_on and depends_on != ".":
                return NXtransformationsAxis(self._handle[depends_on])
        return None

    @cached_property
    def temperature(self) -> pint.Quantity | None:
        """The temperature of the sample."""
        if temperature := self._handle.get("temperature"):
            return temperature[()] * units(temperature)
        return None

    @cached_property
    def transformations(self) -> list[NXtransformations]:
        """This is the recommended location for sample goniometer and other related axes.

        This is a requirement to describe for any scan experiment. The reason it is
        optional is mainly to accommodate XFEL single shot exposures.

        Use of the depends_on field and the NXtransformations group is strongly
        recommended. As noted above this should be an absolute requirement to have for
        any scan experiment.

        The reason it is optional is mainly to accommodate XFEL single shot exposures.
        """
        return [
            NXtransformations(transformations)
            for transformations in self._transformations
        ]


class NXinstrument(H5Mapping):
    """Collection of the components of the instrument or beamline.

    Template of instrument descriptions comprising various beamline components. Each
    component will also be a NeXus group defined by its distance from the sample.
    Negative distances represent beamline components that are before the sample while
    positive distances represent components that are after the sample. This device
    allows the unique identification of beamline components in a way that is valid for
    both reactor and pulsed instrumentation.
    """

    def __init__(self, handle):
        super().__init__(handle)

        (
            self._attenuators,
            self._detector_groups,
            self._detectors,
            self._beams,
            self._transformations,
        ) = find_classes(
            handle,
            "NXattenuator",
            "NXdetector_group",
            "NXdetector",
            "NXbeam",
            "NXtransformations",
        )

    @cached_property
    def transformations(self) -> NXtransformations:
        """
        Transformations relating to the diffractometer but not to the sample.

        These might include a rotation to represent a 2θ arm on which a detector is
        mounted, or a translation of the detector.
        """
        return [
            NXtransformations(transformations)
            for transformations in self._transformations
        ]

    @cached_property
    def name(self) -> str:
        """Name of instrument.

        Consistency with the controlled vocabulary beamline naming in
        https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_diffrn_source.pdbx_synchrotron_beamline.html
        and
        https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_diffrn_source.type.html
        is highly recommended.
        """
        return h5str(self._handle["name"][()])

    @cached_property
    def short_name(self) -> str:
        """Short name for instrument, perhaps the acronym."""
        return h5str(self._handle["name"].attrs.get("short_name"))

    @cached_property
    def time_zone(self) -> str | None:
        """ISO 8601 time_zone offset from UTC."""
        return self._handle.get("time_zone")

    @cached_property
    def attenuators(self):
        return self._attenuators

    @cached_property
    def detector_groups(self) -> list[NXdetector_group]:
        """Optional logical grouping of detectors."""
        return [NXdetector_group(group) for group in self._detector_groups]

    @cached_property
    def detectors(self) -> list[NXdetector]:
        """A detector, detector bank, or multidetector.

        Normally the detector group will have the name detector. However, in the case of
        multiple detectors, each detector needs a uniquely named NXdetector.
        """
        return [NXdetector(detector) for detector in self._detectors]

    @cached_property
    def beams(self) -> list[NXbeam]:
        """Properties of the neutron or X-ray beam at a given location."""
        return [NXbeam(beam) for beam in self._beams]


class NXdetector_group(H5Mapping):
    """Optional logical grouping of detectors.

    Each detector is represented as an NXdetector with its own detector data array. Each
    detector data array may be further decomposed into array sections by use of
    NXdetector_module groups. Detectors can be grouped logically together using
    NXdetector_group. Groups can be further grouped hierarchically in a single
    NXdetector_group (for example, if there are multiple detectors at an endstation or
    multiple endstations at a facility). Alternatively, multiple NXdetector_groups can
    be provided.

    The groups are defined hierarchically, with names given in the group_names field,
    unique identifying indices given in the field group_index, and the level in the
    hierarchy given in the group_parent field. For example if an x-ray detector group,
    DET, consists of four detectors in a rectangular array:

        DTL    DTR
        DLL    DLR

    We could have:

        group_names: ["DET", "DTL", "DTR", "DLL", "DLR"]
        group_index: [1, 2, 3, 4, 5]
        group_parent:  [-1, 1, 1, 1, 1]
    """

    @cached_property
    def group_names(self) -> np.ndarray:
        """
        An array of the names of the detectors or the names of hierarchical groupings of
        detectors.
        """
        return self._handle["group_names"].asstr()[()]

    @cached_property
    def group_index(self) -> NXIntT:
        """An array of unique identifiers for detectors or groupings of detectors.

        Each ID is a unique ID for the corresponding detector or group named in the
        field group_names. The IDs are positive integers starting with 1.
        """
        return self._handle["group_index"][()]

    @cached_property
    def group_parent(self) -> NXIntT:
        """
        An array of the hierarchical levels of the parents of detectors or groupings of
        detectors.

        A top-level grouping has parent level -1.
        """
        return self._handle["group_parent"][()]


class NXdetector(H5Mapping):
    """
    A detector, detector bank, or multidetector.

    Normally the detector group will have the name detector. However, in the case of
    multiple detectors, each detector needs a uniquely named NXdetector.
    """

    def __init__(self, handle):
        super().__init__(handle)
        (self._modules,) = find_classes(handle, "NXdetector_module")

    @cached_property
    def depends_on(self) -> NXtransformationsAxis | None:
        """The axis on which the detector position depends.

        NeXus path to the detector positioner axis that most directly supports the
        detector. In the case of a single-module detector, the detector axis chain may
        start here.
        """
        if "depends_on" in self._handle:
            return NXtransformationsAxis(self._handle[self._handle["depends_on"][()]])
        return None

    @cached_property
    def data(self) -> NXNumber | None:
        """The raw data array for this detector.

        For a dimension-2 detector, the rank of the data array will be 3. For a
        dimension-3 detector, the rank of the data array will be 4. This allows for the
        introduction of the frame number as the first index.
        """
        if "data" in self._handle:
            return self._handle["data"][()]
        return None

    @cached_property
    def description(self) -> str | None:
        """name/manufacturer/model/etc. information."""
        if "description" in self._handle:
            return h5str(np.squeeze(self._handle["description"])[()])
        return None

    @cached_property
    def distance(self) -> pint.Quantity | None:
        """Distance from the sample to the beam center.

        Normally this value is for guidance only, the proper geometry can be found
        following the depends_on axis chain, But in appropriate cases where the
        dectector distance to the sample is observable independent of the axis chain,
        that may take precedence over the axis chain calculation.
        """
        if distance := self._handle.get("distance"):
            return np.squeeze(distance[()] * units(distance))
        return None

    @cached_property
    def distance_derived(self) -> bool | None:
        """Boolean to indicate if the distance is a derived, rather than a primary
        observation.

        If distance_derived true or is not specified, the distance is assumed to be
        derived from detector axis specifications.
        """
        if "distance_derived" in self._handle:
            return bool(self._handle["distance_derived"][()])
        return None

    @cached_property
    def count_time(self) -> pint.Quantity | None:
        """Elapsed actual counting time."""
        if count_time := self._handle.get("count_time"):
            return np.squeeze(count_time[()] * units(count_time, default="seconds"))
        return None

    @cached_property
    def beam_center_x(self) -> pint.Quantity | None:
        """This is the x position where the direct beam would hit the detector.

        This is a length and can be outside of the actual detector. The length can be in
        physical units or pixels as documented by the units attribute. Normally, this
        should be derived from the axis chain, but the direct specification may take
        precedence if it is not a derived quantity.
        """
        if beam_centre_x := self._handle.get("beam_center_x"):
            return np.squeeze(beam_centre_x[()] * units(beam_centre_x, "pixels"))
        return None

    @cached_property
    def beam_center_y(self) -> pint.Quantity | None:
        """This is the y position where the direct beam would hit the detector.

        This is a length and can be outside of the actual detector. The length can be in
        physical units or pixels as documented by the units attribute. Normally, this
        should be derived from the axis chain, but the direct specification may take
        precedence if it is not a derived quantity.
        """
        if beam_centre_y := self._handle.get("beam_center_y"):
            return np.squeeze(beam_centre_y[()] * units(beam_centre_y, "pixels"))
        return None

    @cached_property
    def pixel_mask_applied(self) -> bool | None:
        """
        True when the pixel mask correction has been applied in the electronics, false
        otherwise (optional).
        """
        if "pixel_mask_applied" in self._handle:
            return bool(self._handle["pixel_mask_applied"][()])
        return None

    @cached_property
    def pixel_mask(self) -> NXIntT | None:
        """The 32-bit pixel mask for the detector.

        Can be either one mask for the whole dataset (i.e. an array with indices i, j)
        or each frame can have its own mask (in which case it would be an array with
        indices nP, i, j).

        Contains a bit field for each pixel to signal dead, blind, high or otherwise
        unwanted or undesirable pixels. They have the following meaning:

          - bit 0: gap (pixel with no sensor)
          - bit 1: dead
          - bit 2: under-responding
          - bit 3: over-responding
          - bit 4: noisy
          - bit 5: -undefined-
          - bit 6: pixel is part of a cluster of problematic pixels (bit set in addition
                   to others)
          - bit 7: -undefined-
          - bit 8: user defined mask (e.g. around beamstop)
          - bits 9-30: -undefined-
          - bit 31: virtual pixel (corner pixel with interpolated value)

        Normal data analysis software would not take pixels into account when a bit in
        (mask & 0x0000FFFF) is set. Tag bit in the upper two bytes would indicate
        special pixel properties that normally would not be a sole reason to reject the
        intensity value (unless lower bits are set.

        If the full bit depths is not required, providing a mask with fewer bits is
        permissible.

        If needed, additional pixel masks can be specified by including additional
        entries named pixel_mask_N, where N is an integer. For example, a general bad
        pixel mask could be specified in pixel_mask that indicates noisy and dead
        pixels, and an additional pixel mask from experiment-specific shadowing could be
        specified in pixel_mask_2. The cumulative mask is the bitwise OR of pixel_mask
        and any pixel_mask_N entries.

        If provided, it is recommended that it be compressed.
        """
        if "pixel_mask" in self._handle:
            return self._handle["pixel_mask"][()]
        return None

    @cached_property
    def bit_depth_readout(self) -> int | None:
        """How many bits the electronics record per pixel (recommended)."""
        if "bit_depth_readout" in self._handle:
            return int(self._handle["bit_depth_readout"][()])
        return None

    @cached_property
    def sensor_material(self) -> str:
        """The name of the material a detector sensor is constructed from.

        At times, radiation is not directly sensed by the detector. Rather, the detector
        might sense the output from some converter like a scintillator. This is the name
        of this converter material."""
        return h5str(np.squeeze(self._handle["sensor_material"])[()])

    @cached_property
    def sensor_thickness(self) -> pint.Quantity:
        thickness = self._handle["sensor_thickness"]
        return np.squeeze(thickness)[()] * units(thickness)

    @cached_property
    def underload_value(self) -> int | None:
        """The lowest value at which pixels for this detector would be reasonably be measured.

        For example, given a saturation_value and an underload_value, the valid pixels
        are those less than or equal to the saturation_value and greater than or equal
        to the underload_value.
        """
        if "underload_value" in self._handle:
            return int(self._handle["underload_value"][()])
        return None

    @cached_property
    def saturation_value(self) -> int | None:
        """The value at which the detector goes into saturation.

        Data above this value is known to be invalid.

        For example, given a saturation_value and an underload_value, the valid pixels
        are those less than or equal to the saturation_value and greater than or equal
        to the underload_value.
        """
        if "saturation_value" in self._handle:
            try:
                return int(self._handle["saturation_value"][()])
            except TypeError as e:
                logger.warning(f"Error extracting {self.path}/saturation_value: {e}")
        return None

    @cached_property
    def modules(self) -> list[NXdetector_module]:
        """The list of NXdetector_modules comprising this NXdetector."""
        return [NXdetector_module(module) for module in self._modules]

    @cached_property
    def type(self) -> str | None:
        """Description of type such as scintillator, ccd, pixel, image plate, CMOS, …"""
        if "type" in self._handle:
            return h5str(np.squeeze(self._handle["type"])[()])
        return None

    @cached_property
    def frame_time(self) -> pint.Quantity | None:
        """This is time for each frame. This is exposure_time + readout time."""
        if frame_time := self._handle.get("frame_time"):
            return np.squeeze(frame_time[()] * units(frame_time))
        return None

    @cached_property
    def serial_number(self) -> str | None:
        """Serial number for the detector."""
        if "serial_number" in self._handle:
            return h5str(np.squeeze(self._handle["serial_number"])[()])
        return None


class NXdetector_module(H5Mapping):
    """Representation of the NXdetector_module class.

    Many detectors consist of multiple smaller modules that are operated in sync and
    store their data in a common dataset. To allow consistent parsing of the
    experimental geometry, this application definiton requires all detectors to define a
    detector module, even if there is only one.

    This group specifies the hyperslab of data in the data array associated with the
    detector that contains the data for this module. If the module is associated with a
    full data array, rather than with a hyperslab within a larger array, then a single
    module should be defined, spanning the entire array.
    """

    @cached_property
    def data_origin(self) -> np.ndarray:
        """The offset of this module into the raw data array.

        A dimension-2 or dimension-3 field which gives the indices of the origin of the
        hyperslab of data for this module in the main area detector image in the parent
        NXdetector module.

        The data_origin is 0-based.

        The frame number dimension (nP) is omitted. Thus the data_origin field for a
        dimension-2 dataset with indices (nP, i, j) will be an array with indices
        (i, j), and for a dimension-3 dataset with indices (nP, i, j, k) will be an
        array with indices (i, j, k).

        The order of indices (i, j or i, j, k) is slow to fast.
        """
        origin = self._handle["data_origin"][()]
        assert not isinstance(origin, int)
        return origin

    @cached_property
    def data_size(self) -> np.ndarray:
        """Two or three values for the size of the module in pixels in each direction.

        Dimensionality and order of indices is the same as for data_origin.
        """
        size = self._handle["data_size"][()]
        # Validate that we aren't the int part of NXIntT
        assert not isinstance(size, int)
        return size

    @cached_property
    def data_stride(self) -> NXIntT | None:
        """Two or three values for the stride of the module in pixels in each direction.

        By default the stride is [1,1] or [1,1,1], and this is the most likely case.
        This optional field is included for completeness.
        """
        if "data_stride" in self._handle:
            return self._handle["data_stride"][()]
        return None

    @cached_property
    def module_offset(self) -> NXtransformationsAxis | None:
        """Offset of the module in regards to the origin of the detector in an arbitrary
        direction.
        """
        if "module_offset" in self._handle:
            return NXtransformationsAxis(self._handle["module_offset"])
        return None

    @cached_property
    def fast_pixel_direction(self) -> NXtransformationsAxis:
        """Values along the direction of fastest varying pixel direction.

        The direction itself is given through the vector attribute.
        """
        return NXtransformationsAxis(self._handle["fast_pixel_direction"])

    @cached_property
    def slow_pixel_direction(self) -> NXtransformationsAxis:
        """Values along the direction of slowest varying pixel direction.

        The direction itself is given through the vector attribute.
        """
        return NXtransformationsAxis(self._handle["slow_pixel_direction"])


class NXsource(H5Mapping):
    """The neutron or x-ray storage ring/facility."""

    @cached_property
    def name(self) -> str:
        """Name of source.

        Consistency with the naming in
        https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_diffrn_source.pdbx_synchrotron_site.html
        controlled vocabulary is highly recommended.
        """
        return h5str(self._handle["name"][()])

    @cached_property
    def short_name(self) -> str | None:
        """Short name for source, perhaps the acronym"""
        return h5str(self._handle["name"].attrs.get("short_name"))


class NXbeam(H5Mapping):
    """Properties of the neutron or X-ray beam at a given location.

    It will be referenced by beamline component groups within the NXinstrument group or
    by the NXsample group. Note that variables such as the incident energy could be
    scalar values or arrays. This group is especially valuable in storing the results of
    instrument simulations in which it is useful to specify the beam profile, time
    distribution etc. at each beamline component. Otherwise, its most likely use is in
    the NXsample group in which it defines the results of the neutron scattering by the
    sample, e.g., energy transfer, polarizations.
    """

    @cached_property
    def incident_wavelength(self) -> pint.Quantity:
        """Wavelength on entering beamline component.

        In the case of a monchromatic beam this is the scalar wavelength.

        Several other use cases are permitted, depending on the presence or absence of
        other incident_wavelength_X fields.

        In the case of a polychromatic beam this is an array of length m of wavelengths,
        with the relative weights in incident_wavelength_weights.

        In the case of a monochromatic beam that varies shot- to-shot, this is an array
        of wavelengths, one for each recorded shot. Here, incident_wavelength_weights
        and incident_wavelength_spread are not set.

        In the case of a polychromatic beam that varies shot-to- shot, this is an array
        of length m with the relative weights in incident_wavelength_weights as a 2D
        array.

        In the case of a polychromatic beam that varies shot-to- shot and where the
        channels also vary, this is a 2D array of dimensions nP by m (slow to fast) with
        the relative weights in incident_wavelength_weights as a 2D array.

        Note, variants are a good way to represent several of these use cases in a
        single dataset, e.g. if a calibrated, single-value wavelength value is available
        along with the original spectrum from which it was calibrated.
        """
        wavelength = self._handle["incident_wavelength"]
        return wavelength[()] * units(wavelength)

    @cached_property
    def flux(self) -> pint.Quantity | None:
        """Flux density incident on beam plane area in photons per second per unit area.

        In the case of a beam that varies in flux shot-to-shot, this is an array of
        values, one for each recorded shot.
        """
        if flux := self._handle.get("flux"):
            return flux[()] * units(flux)
        return None

    @cached_property
    def total_flux(self) -> pint.Quantity | None:
        """Flux incident on beam plane in photons per second.

        In the case of a beam that varies in total flux shot-to-shot, this is an array
        of values, one for each recorded shot.
        """
        if total_flux := self._handle.get("total_flux"):
            return total_flux[()] * units(total_flux)
        return None

    @cached_property
    def incident_beam_size(self) -> pint.Quantity | None:
        """Two-element array of FWHM (if Gaussian or Airy function) or diameters
        (if top hat) or widths (if rectangular) of the beam in the order x, y.
        """
        if beam_size := self._handle.get("incident_beam_size"):
            return beam_size[()] * units(beam_size)
        return None

    @cached_property
    def profile(self) -> str | None:
        """The beam profile, Gaussian, Airy function, top-hat or rectangular.

        The profile is given in the plane of incidence of the beam on the sample.

        Any of these values: Gaussian | Airy | top-hat | rectangular
        """
        if "profile" in self._handle:
            return h5str(self._handle["profile"][()])
        return None


@dataclasses.dataclass(frozen=True)
class DependencyChain(Sequence[NXtransformationsAxis]):
    transformations: list[NXtransformationsAxis]

    def __iter__(self) -> Iterator[NXtransformationsAxis]:
        return iter(self.transformations)

    @overload
    def __getitem__(self, idx: int) -> NXtransformationsAxis:
        ...

    @overload
    def __getitem__(self, idx: slice) -> Sequence[NXtransformationsAxis]:
        ...

    def __getitem__(self, idx):
        return self.transformations[idx]

    def __len__(self) -> int:
        return len(self.transformations)

    def __str__(self):
        string = []
        for t in self.transformations:
            depends_on = t.depends_on.path if t.depends_on else "."
            string.extend(
                [
                    f"{t.path} = {t[()]:g}",
                    f"  @transformation_type = {t.transformation_type}",
                    f"  @vector = {t.vector}",
                    f"  @offset = {t.offset}",
                    f"  @depends_on = {depends_on}",
                ]
            )
        return "\n".join(string)


def get_dependency_chain(
    transformation: NXtransformationsAxis,
) -> DependencyChain:
    """Return the dependency chain for a given NXtransformationsAxis.

    Follow the `depends_on` tree for a given NXtransformationsAxis and construct the
    resulting list of NXtransformationsAxis.
    """
    transformations = []
    transform: NXtransformationsAxis | None = transformation
    while transform is not None:
        transformations.append(transform)
        transform = transform.depends_on
    return DependencyChain(transformations)


def get_cumulative_transformation(
    dependency_chain: DependencyChain | Sequence[NXtransformationsAxis],
) -> np.ndarray:
    """Compute the cumulative transformation for a given dependency chain"""
    return reduce(operator.__matmul__, reversed([t.matrix for t in dependency_chain]))


Axes = namedtuple("Axes", ["axes", "angles", "names", "is_scan_axis"])


def get_rotation_axes(dependency_chain: DependencyChain) -> Axes:
    axes = []
    angles = []
    axis_names = []
    is_scan_axis = []

    for transformation in dependency_chain:
        if transformation.transformation_type != "rotation":
            continue
        values = transformation[()].to("degrees").magnitude
        is_scan = len(values) > 1 and not np.all(values == values[0])
        axes.append(transformation.vector)
        angles.append(values[0])
        assert transformation.path
        axis_names.append(transformation.path.split("/")[-1])
        is_scan_axis.append(is_scan)

    return Axes(
        np.array(axes), np.array(angles), np.array(axis_names), np.array(is_scan_axis)
    )
