from __future__ import annotations

import xml
from typing import Dict, Generic, List, Tuple, TypeVar

import numpy as np
from scipy.constants import Planck, m_n
from scipy.interpolate import interp1d

import cctbx.array_family.flex as flex

Shape = TypeVar("Shape")
DType = TypeVar("DType")


class Array(np.ndarray, Generic[Shape, DType]):
    pass


vec3float = Array["3", float]
vec2float = Array["2", float]
vec3int = Array["3", int]
vec2int = Array["2", int]


def wavelength_from_tof(
    distance: float | flex.double, tof: float | flex.double
) -> float | flex.double:
    """
    distance (m)
    tof (s)

    return (A)
    """
    return ((Planck * tof) / (m_n * distance)) * 10**10


def tof_from_wavelength(
    distance: float | flex.double, wavelength: float | flex.double
) -> float | flex.double:
    """
    wavelength (A)
    return (s)
    """

    wavelength = wavelength * 10**-10

    return (wavelength * m_n * distance) / Planck


def frame_to_tof_interpolator(frames: list[float], tof: list[float]) -> interp1d:
    """
    ToF can vary nonlinearly with frame number.
    A cubic spline is used to better account for these cases.
    """
    assert min(frames) >= 0
    assert min(tof) >= 0
    assert len(frames) == len(tof)
    assert all(i < j for i, j in zip(frames, frames[1:]))
    assert all(i < j for i, j in zip(tof, tof[1:]))
    return interp1d(frames, tof, kind="cubic")


def tof_to_frame_interpolator(tof: list[float], frames: list[float]) -> interp1d:
    """
    ToF can vary nonlinearly with frame number.
    A cubic spline is used to better account for these cases.
    """
    assert min(frames) >= 0
    assert min(tof) >= 0
    assert len(frames) == len(tof)
    assert all(i < j for i, j in zip(frames, frames[1:]))
    assert all(i < j for i, j in zip(tof, tof[1:]))
    return interp1d(tof, frames, kind="cubic")


class InstrumentDefinitionReader:
    """
    Class to obtain instrument information in dxtbx format from an
    Instrument Definition File as described here:
    https://docs.mantidproject.org/nightly/concepts/InstrumentDefinitionFile.html#instrumentdefinitionfile

    Example usage:
    >>> import h5py
    >>> h5_file = h5py.File("path/to/nexus_file")
    >>> xml_string = h5_file["path/to/xml_file"]
    >>> import xml.etree.ElementTree as ET
    >>> xml_file = ET.fromstring(xml_string)
    >>> from dxtx.model.tof_helpers import InstrumentDefinitionReader
    >>> reader = InstrumentDefinitionReader()
    >>> origins, fast_axes, slow_axes = reader.get_dials_detector_geometry(xml_file)
    """

    class Panel:
        """
        Class to store panel properties in spherical coordinates
        """

        def __init__(
            self,
            idx: int,
            centre_origin_in_m: float,
            gam_in_deg: float,
            nu_in_deg: float,
            num_pixels: vec2int,
            pixel_size_in_m: vec2float,
            x_orientation: vec2int = np.array((1, 0)),
            y_orientation: vec2int = np.array((0, 1)),
        ) -> None:
            self.idx = idx
            self.centre_origin_in_m = centre_origin_in_m
            self.gam_in_deg = gam_in_deg
            self.nu_in_deg = nu_in_deg
            self.num_pixels = num_pixels
            self.pixel_size_in_m = pixel_size_in_m
            self.x_orientation = x_orientation
            self.y_orientation = y_orientation

        def panel_size_in_m(self) -> vec2float:
            panel_size = []
            for i in range(len(self.num_pixels)):
                panel_size.append(self.num_pixels[i] * self.pixel_size_in_m[i])
            return np.array(panel_size)

        def orientations_flipped(self) -> bool:
            abs_x = tuple([abs(i) for i in self.x_orientation])
            abs_y = tuple([abs(i) for i in self.y_orientation])
            return abs_x == (0, 1) and abs_y == (1, 0)

        def orientation_direction_flipped(self) -> bool:
            return sum(self.x_orientation) == sum(self.y_orientation) == -1

        def __repr__(self) -> None:
            return f"idx: {self.idx} \n \
            centre origin (m): {self.centre_origin_in_m} \n \
            gam (deg): {self.gam_in_deg} \n \
            nu (deg): {self.nu_in_deg} \n \
            num pixels: {self.num_pixels} \n \
            pixel size (m): {self.pixel_size_in_m} \n \
            x_orientation: {self.x_orientation} \n \
            y_orientation: {self.y_orientation}"

    def is_panel(self, tree_component: xml.etree.ElementTree.Element) -> bool:
        if "type" not in tree_component.attrib:
            return False
        return (
            "panel" in tree_component.attrib["type"]
            and "location" in tree_component[0].tag
        )

    def get_rotation_vals(
        self, rot: xml.etree.ElementTree.Element
    ) -> Tuple[float, vec3float]:
        if "rot" in rot.attrib:
            val = float(rot.attrib["rot"])
        elif "val" in rot.attrib:
            val = float(rot.attrib["val"])
        else:
            return None
        try:
            x = int(rot.attrib["axis-x"])
            y = int(rot.attrib["axis-y"])
            z = int(rot.attrib["axis-z"])
        except KeyError:  # if no axes given axis-z is assumed
            x = 0
            y = 0
            z = 1
        return (val, (x, y, z))

    def get_rotations(
        self,
        line: xml.etree.ElementTree.Element,
        rotations: List[Tuple[float, vec3float]],
    ) -> List[Tuple[float, vec3float]]:
        rotation_vals = self.get_rotation_vals(line)
        if rotation_vals is not None:
            rotations.append(self.get_rotation_vals(line))
        try:
            return self.get_rotations(line[0], rotations=rotations)
        except IndexError:
            return rotations

    def x_project_v(self, x: vec3float, v: vec3float) -> Dict[str, vec3float]:
        """Project x onto v, returning parallel and perpendicular components
        >> d = xProject(x, v)
        >> np.allclose(d['par'] + d['perp'], x)
        True
        """
        par = self.x_par_v(x, v)
        perp = x - par
        return {"par": par, "perp": perp}

    def rotate_about(self, a: vec3float, b: vec3float, theta: float) -> vec3float:
        if np.isclose(abs(np.dot(a, b)), 1):
            return a
        """Rotate vector a about vector b by theta radians."""

        proj = self.x_project_v(a, b)
        w = np.cross(b, proj["perp"])
        unit_w = w / np.linalg.norm(w)
        return (
            proj["par"]
            + proj["perp"] * np.cos(theta)
            + np.linalg.norm(proj["perp"]) * unit_w * np.sin(theta)
        )

    def x_par_v(self, x: vec3float, v: vec3float) -> vec3float:
        """Project x onto v. Result will be parallel to v."""
        return np.dot(x, v) / np.dot(v, v) * v

    def x_perp_v(self, x: vec3float, v: vec3float) -> vec3float:
        """Component of x orthogonal to v. Result is perpendicular to v."""
        return x - self.x_par_v(x, v)

    def get_panel_axes(
        self,
        xml_file: xml.etree.ElementTree.Element,
        init_fast_axis: vec3float,
        init_slow_axis: vec3float,
    ) -> Tuple[vec3float, vec3float, vec3float]:
        def get_rot_axis(
            raw_rot_axis: vec3float,
            fast_axis: vec3float,
            slow_axis: vec3float,
            z_axis: vec3float,
            init_fast_axis: vec3float,
            init_slow_axis: vec3float,
            init_z_axis: vec3float,
        ) -> vec3float:
            if abs(np.dot(raw_rot_axis, init_fast_axis)) == 1:
                return fast_axis
            if abs(np.dot(raw_rot_axis, init_slow_axis)) == 1:
                return slow_axis
            if abs(np.dot(raw_rot_axis, init_z_axis)) == 1:
                return z_axis
            raise NotImplementedError

        def get_axes(
            rotations: List[Tuple[float, vec3float]],
            init_fast_axis: vec3float,
            init_slow_axis: vec3float,
            init_z_axis: vec3float,
        ) -> Tuple[Tuple[vec3float], Tuple[vec3float]]:
            fast_axis = init_fast_axis
            slow_axis = init_slow_axis
            z_axis = init_z_axis
            for rot in rotations:
                raw_rot_axis = np.array(rot[1])
                rot_val = np.deg2rad(rot[0])
                rot_axis = get_rot_axis(
                    raw_rot_axis,
                    fast_axis,
                    slow_axis,
                    z_axis,
                    init_fast_axis,
                    init_slow_axis,
                    init_z_axis,
                )

                fast_axis = self.rotate_about(fast_axis, rot_axis, rot_val)
                slow_axis = self.rotate_about(slow_axis, rot_axis, rot_val)
                z_axis = self.rotate_about(z_axis, rot_axis, rot_val)
            return fast_axis, slow_axis, z_axis

        slow_axes = []
        fast_axes = []
        init_z_axis = np.cross(init_fast_axis, init_slow_axis)
        for child in xml_file:
            if self.is_panel(child):
                panel = child[0]
                rotations = []
                rotations = self.get_rotations(
                    line=panel,
                    rotations=rotations,
                )
                axes = get_axes(rotations, init_fast_axis, init_slow_axis, init_z_axis)
                fast_axes.append(tuple(axes[0]))
                slow_axes.append(tuple(axes[1]))
        return fast_axes, slow_axes

    def shift_origin_centre_to_top_left(
        self,
        centre_origin_in_m: vec3float,
        fast_axis: vec3float,
        slow_axis: vec3float,
        panel_size_in_m: vec2float,
    ) -> vec3float:
        # Assumes fast and slow axes are unit vectors

        slow_axis = slow_axis * panel_size_in_m[1] * 0.5
        fast_axis = fast_axis * panel_size_in_m[0] * 0.5
        return centre_origin_in_m - slow_axis - fast_axis

    def rotations_to_spherical_coordinates(
        self,
        zeroth_pixel_origin: vec2float,
        rotations: Tuple[Tuple[float, vec3int], ...],
    ) -> vec2float:
        """
        Corrects rotations for zeroth_pixel_origin and returns
        the gam and nu angles in spherical coordinates.
        """

        gam = rotations[0][0]

        try:
            nu = rotations[1][0]
        except IndexError:
            nu = 0

        xstart, ystart = zeroth_pixel_origin
        if np.sign(xstart) == np.sign(ystart):
            if abs(nu) > 0:
                nu *= -1
            else:
                gam *= -1

        elif abs(nu) > 0:
            gam -= 180

        return gam, nu

    def get_panel_types(
        self, xml_file: xml.etree.ElementTree.Element
    ) -> Dict[str, float | int]:
        def is_panel_settings(tree_element: xml.etree.ElementTree.Element) -> bool:
            required_fields = [
                "xstart",
                "ystart",
                "xpixels",
                "ypixels",
                "xstep",
                "ystep",
            ]
            for i in required_fields:
                if i not in tree_element.attrib:
                    return False

            return True

        panel_types = {}

        for child in xml_file:
            if is_panel_settings(child):
                key = child.attrib["name"]
                xstart = float(child.attrib["xstart"])
                ystart = float(child.attrib["ystart"])
                xpixels = int(child.attrib["xpixels"])
                ypixels = int(child.attrib["ypixels"])
                xpixel_size = abs(float(child.attrib["xstep"]))
                ypixel_size = abs(float(child.attrib["ystep"]))
                panel_types[key] = {
                    "xstart": xstart,
                    "ystart": ystart,
                    "xpixels": xpixels,
                    "ypixels": ypixels,
                    "xpixel_size": xpixel_size,
                    "ypixel_size": ypixel_size,
                }

        return panel_types

    def get_panels(self, xml_file: xml.etree.ElementTree.Element) -> Tuple[Panel]:
        def panel_name_to_idx(name: str) -> int:
            return int("".join(filter(str.isdigit, name)))

        def get_cartesian_origin(panel_attrib):
            try:
                x = float(panel_attrib["x"])
                y = float(panel_attrib["y"])
                z = float(panel_attrib["z"])
            except KeyError:  # assume spherical coordinates
                r = float(panel_attrib["r"])
                t = np.deg2rad(float(panel_attrib["t"]))
                p = np.deg2rad(float(panel_attrib["p"]))
                x = r * np.sin(t) * np.cos(p)
                y = r * np.sin(t) * np.sin(p)
                z = r * np.cos(t)
            return x, y, z

        def get_rotations_from_origin(
            panel_attrib: Dict[str, str],
        ) -> Tuple[float, float]:
            return float(panel_attrib["t"]), float(panel_attrib["p"])

        def get_rotation_vals(
            rot: xml.etree.ElementTree.Element,
        ) -> Tuple[float, vec3float]:
            val = float(rot.attrib["val"])
            try:
                x = int(rot.attrib["axis-x"])
                y = int(rot.attrib["axis-y"])
                z = int(rot.attrib["axis-z"])
            except KeyError:  # if no axes given axis-z is assumed
                x = 0
                y = 0
                z = 1
            return (val, (x, y, z))

        def get_rotations(
            line: xml.etree.ElementTree.Element,
            rotations: List[Tuple[float, vec3float]],
        ) -> List[Tuple[float, vec3float]]:
            rotations.append(get_rotation_vals(line))
            try:
                return get_rotations(line[0], rotations=rotations)
            except IndexError:
                return rotations

        def get_panel_orentation(name: str) -> Tuple[vec2int, vec2int]:
            """
            This is a hack for SXD, where the Mantid SXD_Definition.xml
            does not seem to identify SXD panel 1 as upsidedown.
            """

            if name == "bank1":
                return np.array((0, -1)), np.array((-1, 0))
            else:
                return np.array((1, 0)), np.array((0, 1))

        panels = []
        panel_types = self.get_panel_types(xml_file)

        for child in xml_file:
            if self.is_panel(child):
                panel = child[0]
                name = panel.attrib["name"]
                idx = panel_name_to_idx(name)
                x, y, z = get_cartesian_origin(panel.attrib)
                rotation_start = panel[0]
                rotations = []
                rotations = get_rotations(
                    line=rotation_start,
                    rotations=rotations,
                )

                panel_type = child.attrib["type"]
                panel_info = panel_types[panel_type]

                zeroth_pixel_origin = (panel_info["xstart"], panel_info["ystart"])
                gam_in_deg, nu_in_deg = self.rotations_to_spherical_coordinates(
                    zeroth_pixel_origin=zeroth_pixel_origin, rotations=rotations
                )
                try:
                    gam_in_deg, nu_in_deg = get_rotations_from_origin(panel.attrib)
                except KeyError:
                    gam_in_deg = 0
                    nu_in_deg = 0

                num_pixels = (
                    panel_info["xpixels"],
                    panel_info["ypixels"],
                )

                pixel_size_in_m = (panel_info["xpixel_size"], panel_info["ypixel_size"])

                x_or, y_or = get_panel_orentation(name=name)

                panels.append(
                    InstrumentDefinitionReader.Panel(
                        idx=idx,
                        centre_origin_in_m=(x, y, z),
                        gam_in_deg=gam_in_deg,
                        nu_in_deg=nu_in_deg,
                        num_pixels=num_pixels,
                        pixel_size_in_m=pixel_size_in_m,
                        x_orientation=x_or,
                        y_orientation=y_or,
                    )
                )

        return tuple(panels)

    def panel_idx_to_name(self, idx: int) -> str:
        return "bank" + str(idx)

    def get_panel_names(self, xml_file: xml.etree.ElementTree.Element) -> Tuple[str]:
        panels = self.get_panels(xml_file)
        panel_names = [self.panel_idx_to_name(i.idx) for i in panels]
        return panel_names

    def get_num_panels(self, xml_file: xml.etree.ElementTree.Element) -> int:
        return len(self.get_panels(xml_file))

    def get_dials_detector_geometry(
        self, xml_file: xml.etree.ElementTree.Element
    ) -> Tuple[Tuple[vec3float], Tuple[vec3float], Tuple[vec3float]]:
        # Get panel data in spherical coordinates
        panels = self.get_panels(xml_file)
        panel_dict = {self.panel_idx_to_name(i.idx): i for i in panels}

        # Get panel axes
        init_slow_axis = np.array((0, 1, 0))
        init_fast_axis = np.array((1, 0, 0))

        slow_axes, fast_axes = self.get_panel_axes(
            xml_file, init_fast_axis, init_slow_axis
        )

        # Get panel origins
        origins = []
        count = 0
        for name in panel_dict:
            new_panel = panel_dict[name]
            fast_axis = np.array(fast_axes[count])
            slow_axis = np.array(slow_axes[count])
            top_left_origin_in_m = self.shift_origin_centre_to_top_left(
                centre_origin_in_m=new_panel.centre_origin_in_m,
                fast_axis=fast_axis,
                slow_axis=slow_axis,
                panel_size_in_m=new_panel.panel_size_in_m(),
            )
            top_left_origin_in_mm = np.array([i * 1000 for i in top_left_origin_in_m])

            origins.append(tuple(top_left_origin_in_mm))
            count += 1

        origins = tuple(tuple(float(v) for v in vec) for vec in origins)
        fast_axes = tuple(tuple(float(v) for v in vec) for vec in fast_axes)
        slow_axes = tuple(tuple(float(v) for v in vec) for vec in slow_axes)
        return origins, fast_axes, slow_axes
