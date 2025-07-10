from __future__ import annotations

import logging

import h5py
import numpy as np

import cctbx.array_family.flex as flex

import dxtbx_flumpy as flumpy
from dxtbx import IncorrectFormatError
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.model import Detector
from dxtbx.model.beam import BeamFactory, PolychromaticBeam, Probe
from dxtbx.model.goniometer import Goniometer, GoniometerFactory
from dxtbx.model.scan import Scan, ScanFactory

logger = logging.getLogger(__name__)


class FormatESSNMX(FormatHDF5):
    """
    Class to read files from NMX
    https://europeanspallationsource.se/instruments/nmx
    preprocessed files in scipp to obtain binned data
    """

    def __init__(self, image_file, **kwargs) -> None:
        if not FormatESSNMX.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self._nxs_file = h5py.File(image_file, "r")
        self._raw_data = None

    @staticmethod
    def understand(image_file: str) -> bool:
        try:
            return FormatESSNMX.is_nmx_file(image_file)
        except (OSError, KeyError):
            return False

    @staticmethod
    def is_nmx_file(image_file: str) -> bool:
        def get_name(image_file):
            try:
                with h5py.File(image_file, "r") as handle:
                    return handle["NMX_data"].attrs["name"]
            except (OSError, KeyError, AttributeError):
                return ""

        return get_name(image_file) == "NMX"

    def get_instrument_name(self) -> str:
        return "NMX"

    def get_experiment_description(self) -> str:
        return "Simulated data"

    def _load_raw_data(self) -> None:
        raw_data = []
        image_size = self._get_image_size()
        total_pixels = image_size[0] * image_size[1]
        num_images = self.get_num_images()
        for i in range(self._get_num_panels()):
            spectra = self._nxs_file["NMX_data"]["detector_1"]["counts"][
                0, total_pixels * i : total_pixels * (i + 1), :
            ]
            spectra = np.reshape(spectra, (image_size[0], image_size[1], num_images))
            raw_data.append(flumpy.from_numpy(spectra))

        self._raw_data = tuple(raw_data)

    def get_raw_data(
        self, index: int, use_loaded_data: bool = False
    ) -> tuple[flex.int]:
        raw_data = []
        image_size = self._get_image_size()
        total_pixels = image_size[0] * image_size[1]

        if use_loaded_data:
            if self._raw_data is None:
                self._load_raw_data()
            for panel in self._raw_data:
                data = panel[:, :, index : index + 1]
                data.reshape(flex.grid(panel.all()[0], panel.all()[1]))
                data.matrix_transpose_in_place()
                raw_data.append(data)

        else:
            for i in range(self._get_num_panels()):
                spectra = self._nxs_file["NMX_data"]["detector_1"]["counts"][
                    0, total_pixels * i : total_pixels * (i + 1), index : index + 1
                ]
                spectra = np.reshape(spectra, image_size)
                raw_data.append(flumpy.from_numpy(np.ascontiguousarray(spectra)))

        return tuple(raw_data)

    def _get_time_channel_bins(self) -> list[float]:
        # (usec)
        return self._nxs_file["NMX_data"]["detector_1"]["t_bin"][:] * 10**6

    def _get_time_of_flight(self) -> list[float]:
        # (usec)
        bins = self._get_time_channel_bins()
        return [float((bins[i] + bins[i + 1]) * 0.5) for i in range(len(bins) - 1)]

    def get_num_images(self) -> int:
        return len(self._get_time_of_flight())

    def get_detector(self, index: int = None) -> Detector:
        num_panels = self._get_num_panels()
        panel_names = self._get_panel_names()
        panel_type = self._get_panel_type()
        image_size = self._get_image_size()
        trusted_range = self._get_panel_trusted_range()
        pixel_size = self._get_pixel_size()
        fast_axes = self._get_panel_fast_axes()
        slow_axes = self._get_panel_slow_axes()
        panel_origins = self._get_panel_origins()
        gain = self._get_panel_gain()
        panel_projections = self._get_panel_projections_2d()
        detector = Detector()
        root = detector.hierarchy()

        for i in range(num_panels):
            panel = root.add_panel()
            panel.set_type(panel_type)
            panel.set_name(panel_names[i])
            panel.set_image_size(image_size)
            panel.set_trusted_range(trusted_range)
            panel.set_pixel_size(pixel_size)
            panel.set_local_frame(fast_axes[i], slow_axes[i], panel_origins[i])
            panel.set_gain(gain)
            r, t = panel_projections[i]
            r = tuple(map(int, r))
            t = tuple(map(int, t))
            panel.set_projection_2d(r, t)

        return detector

    def _get_num_panels(self) -> int:
        return self._nxs_file["NMX_data/instrument"].attrs["nr_detector"]

    def _get_panel_names(self) -> list[str]:
        return [
            "%02d" % (i + 1)
            for i in range(self._nxs_file["NMX_data/instrument"].attrs["nr_detector"])
        ]

    def _get_panel_type(self) -> str:
        return "SENSOR_PAD"

    def _get_image_size(self) -> tuple[int, int]:
        # (px)
        return (1280, 1280)

    def _get_panel_trusted_range(self) -> tuple[int, int]:
        # 4 * 1280**2 plus buffer
        return (-1, 7000000)

    def _get_pixel_size(self) -> tuple[float, float]:
        # (mm)
        return (0.4, 0.4)

    def _get_panel_fast_axes(self) -> tuple[tuple[float, float, float]]:
        return ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, -1.0))

    def _get_panel_slow_axes(self) -> tuple[tuple[float, float, float]]:
        return ((0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0))

    def _get_panel_origins(self) -> tuple[tuple[float, float, float]]:
        # (mm)
        return ((-250, -250.0, -292.0), (290, -250.0, -250), (-290, -250.0, 250.0))

    def _get_panel_projections_2d(self) -> dict[int : tuple[tuple, tuple]]:
        p_w, p_h = self._get_image_size()
        p_w += 10
        p_h += 10
        panel_pos = {
            0: ((-1, 0, 0, -1), (p_h, 0)),
            1: ((-1, 0, 0, -1), (p_h, p_w)),
            2: ((-1, 0, 0, -1), (p_h, -p_w)),
        }

        return panel_pos

    def get_beam(self, index: int = None) -> PolychromaticBeam:
        direction = self._get_sample_to_source_direction()
        distance = self._get_sample_to_source_distance()
        wavelength_range = self._get_wavelength_range()
        return BeamFactory.make_polychromatic_beam(
            direction=direction,
            sample_to_source_distance=distance,
            probe=Probe.neutron,
            wavelength_range=wavelength_range,
        )

    def _get_sample_to_source_direction(self) -> tuple[float, float, float]:
        return (0, 0, 1)

    def _get_wavelength_range(self) -> tuple[float, float]:
        # (A)
        return (1.8, 2.55)

    def _get_sample_to_source_distance(self) -> float:
        try:
            dist = abs(self._nxs_file["NMX_data/NXsource/distance"][...]) * 1000
            return dist
        except (KeyError, ValueError):
            logger.warning("sample to moderator_distance not found, using dummy value")
            return 157406

    def _get_panel_gain(self) -> float:
        return 1.0

    def get_goniometer_phi_angle(self) -> float:
        return self.get_goniometer_orientations()[0]

    def get_goniometer(self, index: int = None) -> Goniometer:
        rotation_axis = (0.0, 1.0, 0.0)
        fixed_rotation = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        goniometer = GoniometerFactory.make_goniometer(rotation_axis, fixed_rotation)
        try:
            angles = self.get_goniometer_orientations()
        except KeyError:
            return goniometer
        axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        for idx, angle in enumerate(angles):
            goniometer.rotate_around_origin(axes[idx], -angle)
        return goniometer

    def get_goniometer_orientations(self) -> tuple[float, float, float]:
        # Angles in deg along x, y, z
        return self._nxs_file["NMX_data/crystal_orientation"][...]

    def get_scan(self, index=None) -> Scan:
        image_range = (1, self.get_num_images())
        properties = {"time_of_flight": tuple(self._get_time_of_flight())}
        return ScanFactory.make_scan_from_properties(
            image_range=image_range, properties=properties
        )

    def get_flattened_data(
        self, image_range: None | tuple = None, scale_data: bool = True
    ) -> tuple[flex.int]:
        """
        Image data summed along the time-of-flight direction
        """

        panel_size = self._get_image_size()
        total_pixels = panel_size[0] * panel_size[1]
        max_val = None
        num_tof_bins = len(self._get_time_channel_bins()) - 1
        raw_data = []
        for panel_idx in range(self._get_num_panels()):
            panel_data = self._nxs_file["NMX_data"]["detector_1"]["counts"][
                0, total_pixels * panel_idx : total_pixels * (panel_idx + 1), :
            ]
            panel_data = np.reshape(
                panel_data, (panel_size[0], panel_size[1], num_tof_bins)
            )
            if image_range is not None:
                assert len(image_range) == 2, (
                    "expected image_range to be only two values"
                )
                assert image_range[0] >= 0 and image_range[0] < image_range[1], (
                    "image_range[0] out of range"
                )
                assert image_range[1] <= num_tof_bins, "image_range[1] out of range"
                panel_data = np.sum(
                    panel_data[:, :, image_range[0] : image_range[1]], axis=2
                ).T

            else:
                panel_data = np.sum(panel_data, axis=2).T
            panel_max_val = np.max(panel_data)
            if max_val is None or max_val < panel_max_val:
                max_val = panel_max_val
            raw_data.append(panel_data)

        if scale_data:
            return tuple([(i / max_val).tolist() for i in raw_data])

        return tuple([i.tolist() for i in raw_data])

    def get_flattened_pixel_data(
        self, panel_idx: int, x: int, y: int
    ) -> tuple[tuple, tuple]:
        time_channels = self._get_time_of_flight()
        panel_size = self._get_image_size()
        height = panel_size[1]
        total_pixels = panel_size[0] * panel_size[1]
        idx = (panel_idx * total_pixels) + panel_idx + x * height + y
        return (
            time_channels,
            tuple(self._nxs_file["NMX_data/detector_1/counts"][0, idx, :].tolist()),
        )

    def get_proton_charge(self) -> float:
        return self._nxs_file["NMX_data"]["proton_charge"][...]
