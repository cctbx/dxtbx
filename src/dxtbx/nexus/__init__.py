import logging
from typing import Optional, Tuple, Union

import numpy as np

from cctbx import eltbx
from scitbx.array_family import flex

import dxtbx.model
from dxtbx.format.nexus import dataset_as_flex

from . import nxmx

logger = logging.getLogger(__name__)


KNOWN_SENSOR_MATERIALS = {
    "Si": "Si",
    "Silicon": "Si",
    "CdTe": "CdTe",
    "GaAs": "GaAs",
}


# Conversion from the McStas coordinate system as used by NeXus to the imgCIF
# coordinate system conventionally used by dxtbx:
#   https://manual.nexusformat.org/design.html#design-coordinatesystem
#   https://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Caxis.html
MCSTAS_TO_IMGCIF = np.diag([-1, 1, -1])


def get_dxtbx_goniometer(nxsample: nxmx.NXsample) -> Optional[dxtbx.model.Goniometer]:
    """Generate a dxtbx goniometer model from an NXsample.

    If the NXsample doesn't have a valid depends_on field, then return None.
    """
    if not nxsample.depends_on:
        return None
    dependency_chain = nxmx.get_dependency_chain(nxsample.depends_on)
    logger.debug("Sample dependency chain: %s", dependency_chain)
    axes = nxmx.get_rotation_axes(dependency_chain)
    if len(axes.axes) == 1:
        return dxtbx.model.GoniometerFactory.make_goniometer(
            MCSTAS_TO_IMGCIF @ axes.axes[0], np.identity(3).flatten()
        )
    else:
        if np.sum(axes.is_scan_axis) == 0:
            # A sequence of still images, choose an arbitrary scan axis
            scan_axis = 0
        else:
            assert np.sum(axes.is_scan_axis) == 1, "only one scan axis is supported"
            scan_axis = int(np.where(axes.is_scan_axis)[0][0])
        return dxtbx.model.GoniometerFactory.make_multi_axis_goniometer(
            flex.vec3_double((MCSTAS_TO_IMGCIF @ axes.axes.T).T),
            flex.double(axes.angles),
            flex.std_string(axes.names),
            scan_axis,
        )


def get_dxtbx_beam(nxbeam: nxmx.NXbeam) -> dxtbx.model.Beam:
    """Generate a dxtbx beam model from an NXbeam."""
    return dxtbx.model.BeamFactory.make_beam(
        sample_to_source=(0, 0, 1),
        wavelength=nxbeam.incident_wavelength.to("angstrom").magnitude,
    )


def get_dxtbx_scan(
    nxsample: nxmx.NXsample, nxdetector: nxmx.NXdetector
) -> Optional[dxtbx.model.Scan]:
    """Generate a dxtbx scan model from an NXsample.

    If the NXsample doesn't have a valid depends_on field, then return None.
    """
    if not nxsample.depends_on:
        return None
    dependency_chain = nxmx.get_dependency_chain(nxsample.depends_on)
    logger.debug("Sample dependency chain: %s", dependency_chain)
    scan_axis = None
    for t in dependency_chain:
        # Find the first varying rotation axis
        if (
            t.transformation_type == "rotation"
            and len(t) > 1
            and not np.all(t[()] == t[0])
        ):
            scan_axis = t
            break

    if scan_axis is None:
        # Fall back on the first varying axis of any type
        for t in dependency_chain:
            if len(t) > 1 and not np.all(t[()] == t[0]):
                scan_axis = t
                break

    if scan_axis is None:
        scan_axis = nxsample.depends_on

    is_rotation = scan_axis.transformation_type == "rotation"
    num_images = len(scan_axis)
    image_range = (1, num_images)

    if is_rotation and num_images > 1:
        oscillation = (
            float(scan_axis[0].to("degree").magnitude),
            float((scan_axis[1] - scan_axis[0]).to("degree").magnitude),
        )
    else:
        oscillation = (
            float(scan_axis[0].to("degree").magnitude) if is_rotation else 0,
            0,
        )

    if nxdetector.frame_time is not None:
        frame_time = nxdetector.frame_time.to("seconds").magnitude
        exposure_times = flex.double(num_images, frame_time)
        epochs = flex.double_range(0, num_images) * frame_time
    else:
        exposure_times = flex.double(num_images, 0)
        epochs = flex.double(num_images, 0)

    return dxtbx.model.Scan(
        image_range,
        oscillation,
        exposure_times,
        epochs,
        batch_offset=0,
        deg=True,
    )


def get_dxtbx_detector(
    nxdetector: nxmx.NXdetector, nxbeam: nxmx.NXbeam
) -> dxtbx.model.Detector:
    """Generate a dxtbx detector model from an NXdetector and NXbeam.

    If the NXdetector contains multiple NXdetector_modules, then a hierarchical detector
    will be returned, else a "flat" detector model with a single panel will be returned
    where there is only a single NXdetector_module.
    """

    detector = dxtbx.model.Detector()

    root: Union[dxtbx.model.Detector, dxtbx.model.DetectorNode]
    if len(nxdetector.modules) > 1:
        root = detector.hierarchy()
    else:
        root = detector

    for module in nxdetector.modules:

        if len(nxdetector.modules) > 1:
            # Set up the detector hierarchy
            if module.fast_pixel_direction.depends_on is not None:
                reversed_dependency_chain = reversed(
                    nxmx.get_dependency_chain(module.fast_pixel_direction.depends_on)
                )
                pg = None
                for i, transformation in enumerate(reversed_dependency_chain):
                    name = transformation.path
                    if pg is None:
                        pg = root
                    pg_names = [child.get_name() for child in pg]
                    if name in pg_names:
                        pg = pg[pg_names.index(name)]
                        continue
                    else:
                        pg = pg.add_group()
                    A = transformation.matrix
                    origin = MCSTAS_TO_IMGCIF @ A[0, :3, 3]
                    fast = (
                        MCSTAS_TO_IMGCIF @ (A @ np.array((-1, 0, 0, 1)))[0, :3] - origin
                    )
                    slow = (
                        MCSTAS_TO_IMGCIF @ (A @ np.array((0, 1, 0, 1)))[0, :3] - origin
                    )
                    pg.set_local_frame(fast, slow, origin)
                    pg.set_name(name)
        else:
            # Use a flat detector model
            pg = root

        if isinstance(pg, dxtbx.model.DetectorNode):
            # Hierarchical detector model
            fast_axis = MCSTAS_TO_IMGCIF @ module.fast_pixel_direction.vector
            slow_axis = MCSTAS_TO_IMGCIF @ module.slow_pixel_direction.vector
            origin = np.array((0.0, 0.0, 0.0))
        else:
            # Flat detector model

            # Apply any rotation components of the dependency chain to the fast axis
            assert module.fast_pixel_direction.depends_on is not None
            fast_axis_depends_on = [
                t
                for t in nxmx.get_dependency_chain(
                    module.fast_pixel_direction.depends_on
                )
                if t.transformation_type == "rotation"
            ]
            if fast_axis_depends_on:
                R = nxmx.get_cumulative_transformation(fast_axis_depends_on)[0, :3, :3]
            else:
                R = np.identity(3)
            fast_axis = MCSTAS_TO_IMGCIF @ R @ module.fast_pixel_direction.vector

            # Apply any rotation components of the dependency chain to the slow axis
            assert module.slow_pixel_direction.depends_on is not None
            slow_axis_depends_on = [
                t
                for t in nxmx.get_dependency_chain(
                    module.slow_pixel_direction.depends_on
                )
                if t.transformation_type == "rotation"
            ]
            if slow_axis_depends_on:
                R = nxmx.get_cumulative_transformation(slow_axis_depends_on)[0, :3, :3]
            else:
                R = np.identity(3)
            slow_axis = MCSTAS_TO_IMGCIF @ R @ module.slow_pixel_direction.vector

            # Apply all components of the dependency chain to the module offset to get the
            # dxtbx panel origin
            dependency_chain = nxmx.get_dependency_chain(
                module.fast_pixel_direction.depends_on
            )
            A = nxmx.get_cumulative_transformation(dependency_chain)
            origin = MCSTAS_TO_IMGCIF @ A[0, :3, 3]

        pixel_size = (
            module.fast_pixel_direction[()].to("mm").magnitude.item(),
            module.slow_pixel_direction[()].to("mm").magnitude.item(),
        )
        # dxtbx requires image size in the order fast, slow - which is the reverse of what
        # is stored in module.data_size
        image_size = tuple(map(int, module.data_size[::-1]))
        underload = (
            float(nxdetector.underload_value)
            if nxdetector.underload_value is not None
            else -0x7FFFFFFF
        )
        overload = (
            float(nxdetector.saturation_value)
            if nxdetector.saturation_value is not None
            else 0x7FFFFFFF
        )
        trusted_range = (underload, overload)

        material = KNOWN_SENSOR_MATERIALS.get(nxdetector.sensor_material)
        if not material:
            raise ValueError(f"Unknown material: {nxdetector.sensor_material}")
        thickness = nxdetector.sensor_thickness.to("mm").magnitude
        table = eltbx.attenuation_coefficient.get_table(material)
        mu = (
            table.mu_at_angstrom(
                nxbeam.incident_wavelength.to("angstrom").magnitude.item()
            )
            / 10.0
        )
        px_mm = dxtbx.model.ParallaxCorrectedPxMmStrategy(mu, thickness)
        name = module.path

        assert pg is not None
        p = pg.add_panel()
        p.set_type("SENSOR_PAD")
        p.set_name(name)
        p.set_local_frame(fast_axis, slow_axis, origin)
        p.set_pixel_size(pixel_size)
        p.set_image_size(image_size)
        p.set_trusted_range(trusted_range)
        p.set_thickness(thickness)
        p.set_material(material)
        p.set_mu(mu)
        p.set_px_mm_strategy(px_mm)

    return detector


def get_detector_module_slices(
    nxdetector: nxmx.NXdetector,
) -> Tuple[Tuple[slice, ...], ...]:
    """Return the slices pointing to the hyperslab of data for each module.

    This will be a tuple of tuples, where each tuple contains the slices corresponding
    to the slow and fast dimensions respectively.
    """
    return tuple(
        tuple(
            slice(int(start), int(start + step), 1)
            for start, step in zip(module.data_origin, module.data_size)
        )
        for module in nxdetector.modules
    )


def get_static_mask(nxdetector: nxmx.NXdetector) -> Tuple[flex.bool, ...]:
    """Return the static mask for an NXdetector.

    This will be a tuple of flex.bool, of length equal to the number of modules. The
    result is intended to be compatible with the get_static_mask() method of dxtbx
    format classes.
    """
    pixel_mask = nxdetector.get("pixel_mask")
    assert pixel_mask and pixel_mask.ndim == 2
    all_slices = get_detector_module_slices(nxdetector)
    return tuple(dataset_as_flex(pixel_mask, slices) == 0 for slices in all_slices)


def get_raw_data(
    nxdata: nxmx.NXdata, nxdetector: nxmx.NXdetector, index: int
) -> Tuple[Union[flex.float, flex.double, flex.int], ...]:
    """Return the raw data for an NXdetector.

    This will be a tuple of flex.float, flex.double or flex.int arrays, of length equal
    to the number of modules. The result is intended to be compatible with the
    get_raw_data() method of dxtbx format classes.
    """
    if nxdata.signal:
        data = nxdata[nxdata.signal]
    else:
        data = list(nxdata.values())[0]
    all_data = []
    for module_slices in get_detector_module_slices(nxdetector):
        slices = [slice(index, index + 1, 1)]
        slices.extend(module_slices)
        data_as_flex = dataset_as_flex(data, tuple(slices))
        # Convert a slice of a 3- or 4-dimension array to a 2D array
        data_as_flex.reshape(flex.grid(data_as_flex.all()[-2:]))
        all_data.append(data_as_flex)
    return tuple(all_data)
