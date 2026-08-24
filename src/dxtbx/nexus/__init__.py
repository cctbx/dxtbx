from __future__ import annotations

import itertools
import logging
import pathlib
from typing import Literal, Optional

import h5py
import numpy as np
import nxmx

import cctbx
from cctbx import eltbx
from scitbx.array_family import flex

import dxtbx.model
from dxtbx import flumpy

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


def get_dxtbx_goniometer(nxsample: nxmx.NXsample) -> dxtbx.model.Goniometer | None:
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


class CachedWavelengthBeamFactory:
    """Defer Beam generation whilst caching the wavelength value"""

    def __init__(self, nxbeam: nxmx.NXbeam):
        self.handle = nxbeam._handle
        self.index = None
        self.model = None
        self.spectrum = None

    def make_beam(self, index: int = 0) -> dxtbx.model.Beam:
        self.read_models(index)
        return self.model

    def make_spectrum(self, index: int = 0) -> dxtbx.model.Spectrum:
        self.read_models(index)
        return self.spectrum

    def get_wavelengths(self) -> list:
        """Read 1-D incident_wavelength array; return per-frame wavelengths in Angstrom."""
        primary_key = "incident_wavelength"
        wavelength = self.handle[primary_key]
        if wavelength.shape in ((), (1,)):
            raise RuntimeError(
                "get_wavelengths requires a 1-D wavelength array; "
                "this file has a scalar incident_wavelength"
            )
        wavelength_units = nxmx.units(wavelength)
        wl_values = wavelength[()]
        return [
            float((w * wavelength_units).to("angstrom").magnitude) for w in wl_values
        ]

    def read_models(self, index: int = 0):
        # Cached model
        if self.model is not None and index == self.index:
            return

        # Get the items from the NXbeam class
        # Note it would be better if there were a general way to read weights
        # and variants from the nxmx classes
        # See https://github.com/cctbx/dxtbx/issues/549
        primary_key = "incident_wavelength"
        wavelength = self.handle[primary_key]
        spectrum_wavelengths = wavelength
        spectrum_weights = self.handle.get(primary_key + "_weights")
        if spectrum_weights is None:
            # Handle deprecation: https://github.com/nexusformat/definitions/issues/837
            spectrum_weights = self.handle.get(primary_key + "_weight")

        # If the wavelength array does not represent spectra, look for spectra
        # in the variant chain
        variant_test = wavelength
        has_variant_spectra = False
        while spectrum_weights is None:
            if "variant" in variant_test.attrs:
                variant_key = variant_test.attrs["variant"]
                variant_wavelengths = self.handle[variant_key]
                variant_weights = self.handle.get(variant_key + "_weights")
                if variant_weights is None:
                    # Handle deprecation: https://github.com/nexusformat/definitions/issues/837
                    variant_weights = self.handle.get(variant_key + "_weight")
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
                wavelength_value = wavelength[()].item()
            else:
                wavelength_value = wavelength[index]
            wavelength_units = nxmx.units(wavelength)
            wavelength_value = float(
                (wavelength_value * wavelength_units).to("angstrom").magnitude
            )
            return wavelength_value

        if spectrum_weights is None:
            # Construct the beam model
            wavelength_value = get_wavelength(wavelength)
            self.model = dxtbx.model.Beam(
                direction=(0, 0, 1), wavelength=wavelength_value
            )
        else:
            self.model = dxtbx.model.Beam()
            self.model.set_direction((0, 0, 1))

            wavelength_units = nxmx.units(spectrum_wavelengths)

            if len(spectrum_wavelengths.shape) > 1:
                spectrum_wavelengths = spectrum_wavelengths[index]
            else:
                spectrum_wavelengths = spectrum_wavelengths[()]
            if len(spectrum_weights.shape) > 1:
                spectrum_weights = spectrum_weights[index]
            else:
                spectrum_weights = spectrum_weights[()]

            spectrum_wavelengths = (
                (spectrum_wavelengths * wavelength_units).to("angstrom").magnitude
            )
            spectrum_energies = cctbx.factor_ev_angstrom / spectrum_wavelengths
            self.spectrum = dxtbx.model.Spectrum(spectrum_energies, spectrum_weights)

            if has_variant_spectra:
                wavelength_value = get_wavelength(wavelength)
                self.model.set_wavelength(wavelength_value)
            else:
                self.model.set_wavelength(self.spectrum.get_weighted_wavelength())


def get_dxtbx_scan(
    nxsample: nxmx.NXsample, nxdetector: nxmx.NXdetector
) -> dxtbx.model.Scan | None:
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

    oscillation = (0, 0)
    if is_rotation:
        start = scan_axis[0].to("degree")
        if scan_axis.end:
            steps = scan_axis.end[()] - scan_axis[()]
        elif num_images > 1:
            steps = np.diff(scan_axis[()])
        else:
            steps = nxmx.ureg.Quantity(0, "degree")

        step = np.median(steps).to("degree")
        if np.any(np.abs(steps - step) > abs(0.1 * step)):
            logger.warning(
                "One or more recorded oscillation widths differ from the median "
                "by more than 10%. The rotation axis of your goniometer may not "
                "have been scanning at a constant speed throughout the data "
                "collection."
            )

        oscillation = tuple(float(f) for f in (start.magnitude, step.magnitude))

    if nxdetector.frame_time is not None:
        frame_time = float(nxdetector.frame_time.to("seconds").magnitude)
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
    nxdetector: nxmx.NXdetector,
    wavelength: float,
    nxdata: nxmx.NXdata | None = None,
) -> dxtbx.model.Detector:
    """Generate a dxtbx detector model from an NXdetector and NXbeam.

    If the NXdetector contains multiple NXdetector_modules, then a hierarchical detector
    will be returned, else a "flat" detector model with a single panel will be returned
    where there is only a single NXdetector_module.
    """

    detector = dxtbx.model.Detector()

    root: dxtbx.model.Detector | dxtbx.model.Panel
    if len(nxdetector.modules) > 1:
        root = detector.hierarchy()
    else:
        root = detector

    for module in nxdetector.modules:
        if len(nxdetector.modules) > 1:
            # Set up the detector hierarchy
            if module.fast_pixel_direction.depends_on is not None:
                reversed_dependency_chain = list(
                    reversed(
                        nxmx.get_dependency_chain(
                            module.fast_pixel_direction.depends_on
                        )
                    )
                )
                pg: dxtbx.model.Detector | dxtbx.model.Panel | None = None

                # Verify that equipment_components in the dependency chain are
                # 1) contiguous and 2) unique
                found = []
                for dependency in reversed_dependency_chain:
                    if dependency.equipment_component:
                        if dependency.equipment_component in found:
                            assert dependency.equipment_component == found[-1]
                        found.append(dependency.equipment_component)

                # Group any transformations together that share the same equipment_component
                # to reduce the number of hierarchy levels

                # Keep transformations without equipment_component set separate by using
                # a different key
                counter = 0

                def equipment_component_key(dependency):
                    if dependency.equipment_component:
                        return dependency.equipment_component  # always a string
                    else:
                        nonlocal counter
                        counter += 1
                        return counter

                for _, transformation_group in itertools.groupby(
                    reversed_dependency_chain, key=equipment_component_key
                ):
                    transformation_group = list(transformation_group)
                    name = transformation_group[-1].path
                    if pg is None:
                        pg = root
                    else:
                        assert isinstance(
                            pg, (dxtbx.model.Detector, dxtbx.model.DetectorNode)
                        )
                        pg_names = [child.get_name() for child in pg]
                        if name in pg_names:
                            pg = pg[
                                pg_names.index(name)
                            ]  # Getitem always returns panel
                            continue
                        else:
                            pg = pg.add_group()
                    A = nxmx.get_cumulative_transformation(transformation_group)
                    origin = MCSTAS_TO_IMGCIF @ A[0, :3, 3]
                    fast = (
                        MCSTAS_TO_IMGCIF @ (A @ np.array((-1, 0, 0, 1)))[0, :3] - origin
                    )
                    slow = (
                        MCSTAS_TO_IMGCIF @ (A @ np.array((0, 1, 0, 1)))[0, :3] - origin
                    )
                    pg.set_local_frame(fast, slow, origin)
                    assert name is not None
                    pg.set_name(name)
                # assert pg is not None
        else:
            # Use a flat detector model
            pg = root

        pixel_size = (
            module.fast_pixel_direction[()].to("mm").magnitude.item(),
            module.slow_pixel_direction[()].to("mm").magnitude.item(),
        )

        if isinstance(pg, dxtbx.model.DetectorNode):
            # Hierarchical detector model
            fast_axis = MCSTAS_TO_IMGCIF @ module.fast_pixel_direction.vector
            slow_axis = MCSTAS_TO_IMGCIF @ module.slow_pixel_direction.vector
            origin = MCSTAS_TO_IMGCIF @ (
                module.fast_pixel_direction.offset.to("mm").magnitude
                if module.fast_pixel_direction.offset is not None
                else (
                    np.array([0.0, 0.0, 0.0])
                    + module.slow_pixel_direction.offset.to("mm").magnitude
                    if module.slow_pixel_direction.offset is not None
                    else np.array([0.0, 0.0, 0.0])
                )
            )
        else:
            # Flat detector model

            # Apply any rotation components of the dependency chain to the fast axis
            # assert module.fast_pixel_direction.depends_on is not None
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
            # assert module.slow_pixel_direction.depends_on is not None
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
            if module.fast_pixel_direction.depends_on is not None:
                top = module.fast_pixel_direction.depends_on
            else:
                # This works around existing EIGER data putting '.' for depends_on
                top = module.module_offset.depends_on
            dependency_chain = nxmx.get_dependency_chain(top)
            A = nxmx.get_cumulative_transformation(dependency_chain)

            origin = MCSTAS_TO_IMGCIF @ (
                (
                    module.fast_pixel_direction.offset.to("mm").magnitude
                    if module.fast_pixel_direction.offset is not None
                    else (
                        np.array([0.0, 0.0, 0.0])
                        + module.slow_pixel_direction.offset.to("mm").magnitude
                        if module.slow_pixel_direction.offset is not None
                        else np.array([0.0, 0.0, 0.0])
                    )
                )
                + A[0, :3, 3]
            )

            if (
                origin[0] == 0
                and origin[1] == 0
                and nxdetector.beam_center_x
                and nxdetector.beam_center_y
            ):
                # fallback on the explicit beam center - this is needed for some older
                # dectris eiger filewriter datasets
                origin -= nxdetector.beam_center_x.magnitude * pixel_size[0] * fast_axis
                origin -= nxdetector.beam_center_y.magnitude * pixel_size[1] * slow_axis

        # dxtbx requires image size in the order fast, slow - which is the reverse of what
        # is stored in module.data_size. Additionally, data_size can have more than 2
        # dimensions, for multi-module detectors. So take the last two dimensions and reverse
        # them.  Examples:
        # [1,2,3]   --> (3, 2)
        # [1,2]     --> (2, 1)
        image_size = (int(module.data_size[-1]), int(module.data_size[-2]))
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
        # The dxtbx trusted_range is inclusive [min-trusted-value, max-trusted-value]
        trusted_range = (underload, overload)

        material = KNOWN_SENSOR_MATERIALS.get(nxdetector.sensor_material)
        if not material:
            raise ValueError(f"Unknown material: {nxdetector.sensor_material}")
        thickness = float(nxdetector.sensor_thickness.to("mm").magnitude)
        table = eltbx.attenuation_coefficient.get_table(material)
        mu = table.mu_at_angstrom(wavelength) / 10.0
        px_mm = dxtbx.model.ParallaxCorrectedPxMmStrategy(mu, thickness)
        name = module.path

        assert name is not None
        assert pg is not None
        assert isinstance(pg, (dxtbx.model.Detector, dxtbx.model.DetectorNode))
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

        if nxdata and nxdata.data_scale_factor and not nxdata.data_scale_factor.shape:
            p.set_gain(1 / nxdata.data_scale_factor)

    return detector


def get_detector_module_slices(
    nxdetector: nxmx.NXdetector,
) -> tuple[tuple[slice, ...], ...]:
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


def get_static_mask(nxdetector: nxmx.NXdetector) -> tuple[flex.bool, ...] | None:
    """Return the static mask for an NXdetector.

    This will be a tuple of flex.bool, of length equal to the number of modules. The
    result is intended to be compatible with the get_static_mask() method of dxtbx
    format classes.
    """
    try:
        pixel_mask = nxdetector.pixel_mask
    except KeyError:
        return None
    if pixel_mask is None or not pixel_mask.size:
        return None
    all_slices = get_detector_module_slices(nxdetector)
    all_mask_slices = []
    for slices in all_slices:
        mask_slice = flumpy.from_numpy(np.ascontiguousarray(pixel_mask[slices])) == 0
        mask_slice.reshape(
            flex.grid(mask_slice.all()[-2:])
        )  # handle 3 or 4 dimension arrays
        all_mask_slices.append(mask_slice)
    return tuple(all_mask_slices)


def _dataset_as_flex(
    data: h5py.Dataset,
    slices: tuple[slice, ...] | None,
    bit_depth: Literal[32] | None = None,
) -> flex.float | flex.double | flex.int:
    """
    Convert an HDF5 dataset to one of the expected flex types.

    Args:
        data: The HDF5 Dataset to convert
        slices:
            The Dataset will be sliced and made contiguous to this shape.
        bit_depth:
            If set to 32, and the dataset in signed integer
            representation requires more than 32-bit, then the data will
            be converted to a signed 32-bit integer. Without this, such
            attempted conversions will raise a TypeError.
    """
    # Make a guaranteed-contiguous copy of the sliced data
    data_np = np.ascontiguousarray(data[slices or ()])
    dtype = data_np.dtype

    # Handle integer conversion. Safe to convert if:
    # - Is signed and <= 4 bytes
    # - Is unsigned and <= 2 bytes
    #
    # Unsafe conversions to 32-bit integer can occur, but only if
    # bit_depth is explicitly set to 32.
    if np.issubdtype(dtype, np.integer):
        if (
            (np.issubdtype(dtype, np.signedinteger) and dtype.itemsize <= 4)
            or (np.issubdtype(dtype, np.unsignedinteger) and dtype.itemsize <= 2)
            or bit_depth == 32
        ):
            data_np = data_np.astype(np.int32, copy=False)
        else:
            raise TypeError(f"Unsupported integer dtype {data_np.dtype}")
    elif np.issubdtype(dtype, np.floating):
        if dtype.itemsize <= 4:
            # Promote anything <= single precision up to single precision
            data_np = data_np.astype(np.float32, copy=False)
        else:
            # Otherwise, everything else becomes double
            data_np = data_np.astype(np.float64, copy=False)
    else:
        # Isn't a recognised integer or floating point type
        raise TypeError(f"Unsupported dtype {data_np.dtype}")

    data_flex = flumpy.from_numpy(data_np)
    return data_flex


def get_raw_data(
    nxdata: nxmx.NXdata,
    nxdetector: nxmx.NXdetector,
    index: int,
    bit_depth: int | None = None,
    module_slices: tuple[tuple[slice, ...], ...] | None = None,
) -> tuple[flex.float | flex.double | flex.int, ...]:
    """Return the raw data for an NXdetector.

    This will be a tuple of flex.float, flex.double or flex.int arrays, of length equal
    to the number of modules. The result is intended to be compatible with the
    get_raw_data() method of dxtbx format classes.
    """
    if nxdata.signal:
        try:
            data = nxdata[nxdata.signal]
        except KeyError:
            logger.warning(f"Key {nxdata.signal} specified by NXdata.signal missing")
            data = list(nxdata.values())[0]
    else:
        data = list(nxdata.values())[0]
    if module_slices is None:
        module_slices = get_detector_module_slices(nxdetector)
    all_data = []
    sliced_outer = data[index]
    for slices in module_slices:
        data_as_flex = _dataset_as_flex(
            sliced_outer, tuple(slices), bit_depth=bit_depth
        )
        data_as_flex.reshape(
            flex.grid(data_as_flex.all()[-2:])
        )  # handle 3 or 4 dimension arrays
        all_data.append(data_as_flex)
    return tuple(all_data)


def _resolve_source(file_name, master_directory: pathlib.Path) -> pathlib.Path:
    """Resolve a VDS source filename, which may be relative to the master file."""
    src = pathlib.Path(file_name)
    if not src.is_absolute():
        src = master_directory / src
    return src


def _vds_source_map(
    data: h5py.Dataset, master_file: h5py.File, master_directory: pathlib.Path
):
    """Frame-ordered [(global_start, global_end, source_path, dataset_name), ...] for
    the source data files backing a virtual dataset.

    Each source must ultimately resolve to an *external* data file, matching how DLS
    NXmx masters are written. Two layouts are supported:

      * the VDS source names the data file directly (``file_name != "."``); or
      * the VDS source is the master file itself (``file_name == "."``) and names an
        in-master dataset that is an :class:`h5py.ExternalLink` to the data file (the
        DLS/GDA two-level layout).

    A ``"."`` source that resolves to a genuine in-master dataset (no external file) is
    unsupported and raises :class:`ValueError`.
    """
    sources = []
    for vmap in data.virtual_sources():
        start, end = vmap.vspace.get_select_bounds()
        file_name = vmap.file_name
        dset_name = vmap.dset_name
        if file_name == ".":
            # Source lives in the master file; on DLS masters it is an external link
            # to the real data file. Resolve one level to find that file.
            link = master_file.get(dset_name, getlink=True)
            if not isinstance(link, h5py.ExternalLink):
                raise ValueError(
                    f"VDS source {dset_name!r} in {master_file.filename!r} does not "
                    "point to an external data file; unsupported layout for live "
                    "availability checking."
                )
            src_path = _resolve_source(link.filename, master_directory)
            dset_name = link.path
        else:
            src_path = _resolve_source(file_name, master_directory)
        sources.append((int(start[0]), int(end[0]), src_path, dset_name))
    return sorted(sources)


def _written_frame_count(dset) -> int:
    """Number of frames actually written to a source dataset.

    Detector writers allocate the source dataset at its *full* block extent up front
    (``shape[0]`` is the planned frame count from creation, before any pixel is
    written), so the dataset extent cannot tell us how much has arrived. What does
    grow incrementally is chunk allocation: with one chunk per frame
    (``chunks == (1, height, width)``) a chunk is only allocated once its frame has
    been written, so the number of allocated chunks is the number of frames written.

    Falls back to ``shape[0]`` for a genuinely growable/resized dataset or if the
    chunk count cannot be queried (e.g. contiguous storage).
    """
    nframes = int(dset.shape[0])
    if dset.chunks is not None and dset.chunks[0] == 1:
        try:
            n = dset.id.get_num_chunks()
        except (RuntimeError, OSError, ValueError, TypeError):
            return nframes
        # One chunk per frame along the first axis, so allocated chunks == frames
        # written; clamp defensively to the declared extent.
        return min(int(n), nframes)
    return nframes


def _frames_written(src_path: pathlib.Path, dset_name: str, method: str) -> int:
    """Number of frames currently readable in a single VDS source data file.

    Returns 0 if the file/dataset cannot be read yet, e.g. it is still locked by
    the SWMR writer, mid-flush, or not present.
    """
    try:
        if method == "file_close":
            # A plain (non-SWMR) open succeeds only once the writer has released the
            # file, i.e. the whole block has finished; a block still being written is
            # locked and raises (caught below -> 0). Detector files are preallocated at
            # full extent, so count written chunks rather than the (static) extent -- a
            # closed-but-empty pre-created file must not read as fully available.
            with h5py.File(src_path, "r") as fh:
                return _written_frame_count(fh[dset_name])
        else:  # "swmr"
            with h5py.File(src_path, "r", swmr=True) as fh:
                dset = fh[dset_name]
                if hasattr(dset, "refresh"):
                    dset.refresh()
                return _written_frame_count(dset)
    except (OSError, KeyError):
        return 0


def get_frame_counts(master_path, method: str = "swmr") -> tuple[int, int]:
    """Return ``(available, total)`` frame counts for an NXmx master.

    During live data collection the master file's virtual dataset (VDS) declares the
    full planned number of frames up front, while the underlying source data files are
    filled incrementally (typically via SWMR). Reading an unwritten region of the VDS
    silently returns the fill value, so this inspects the VDS sources directly to report
    how many frames are genuinely available.

    * ``available``: number of contiguous leading frames actually written to disk.
    * ``total``: the full planned number of frames the master declares.

    ``method`` selects how a source data file's readiness is judged:
      * ``"swmr"`` (default): count the frames written to each source via an SWMR read
        (per-image granularity). Detector datasets are allocated at their full extent up
        front, so this counts allocated per-frame chunks rather than the (static) extent.
      * ``"file_close"``: treat a source as ready only once it can be opened without
        SWMR, i.e. the writer has closed it (per-block granularity).

    A non-virtual dataset (a completed file, or a format without a VDS) is reported as
    fully available (``available == total``). Returns ``(0, 0)`` if the master itself
    cannot be read.
    """
    master_path = pathlib.Path(master_path)
    try:
        fh = h5py.File(master_path, "r", swmr=True)
    except (OSError, ValueError):
        try:
            fh = h5py.File(master_path, "r")
        except OSError:
            return 0, 0
    try:
        nxmx_obj = nxmx.NXmx(fh)
        nxdata = nxmx_obj.entries[0].data[0]
        data = nxdata[nxdata.signal] if nxdata.signal else list(nxdata.values())[0]
        total = int(data.shape[0])
        if not data.is_virtual:
            return total, total
        available = 0
        for gstart, gend, src, dset_name in _vds_source_map(
            data, fh, master_path.parent
        ):
            if gstart > available:
                break  # gap before this source: no longer contiguous from frame 0
            span = min(gend, total - 1) - gstart + 1
            written = _frames_written(src, dset_name, method)
            available = gstart + min(written, span)
            if written < span:
                break  # this block is only partially written; stop here
        return available, total
    finally:
        fh.close()


def get_available_frame_count(master_path, method: str = "swmr") -> int:
    """Number of contiguous leading frames of an NXmx master actually written to disk.

    Thin wrapper over :func:`get_frame_counts` returning only the available count. See
    that function for the full description of ``method`` and the availability model.
    """
    return get_frame_counts(master_path, method)[0]
