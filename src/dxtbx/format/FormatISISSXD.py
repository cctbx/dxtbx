from __future__ import annotations

from typing import List, Tuple

import h5py

import cctbx.array_family.flex as flex

from dxtbx import IncorrectFormatError
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.model import Detector, Goniometer, Scan
from dxtbx.model.beam import BeamFactory, PolychromaticBeam, Probe
from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model.scan import ScanFactory


class FormatISISSXD(FormatHDF5):
    """
    Class to read image files from ISIS SXD
    (https://www.isis.stfc.ac.uk/Pages/sxd.aspx)
    """

    def __init__(self, image_file: str, **kwargs) -> None:
        super().__init__(image_file, **kwargs)
        if not FormatISISSXD.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self._nxs_file = h5py.File(image_file, "r")
        self._raw_data = None

    @staticmethod
    def understand(image_file: str) -> bool:
        try:
            return FormatISISSXD.is_isissxd_file(image_file)
        except (IOError, KeyError):
            return False

    @staticmethod
    def is_isissxd_file(image_file: str) -> bool:
        def get_name(image_file: str) -> str:
            with h5py.File(image_file, "r") as handle:
                if "raw_data_1" in handle:
                    return handle["/raw_data_1/name"][0].decode()
                return ""

        return get_name(image_file) == "SXD"

    def get_experiment_title(self) -> str:
        return self._nxs_file["raw_data_1"]["title"][0].decode()

    def get_experiment_run_number(self) -> str:
        return self._nxs_file["raw_data_1"]["run_number"][0]

    def get_experiment_description(self) -> str:
        title = self.get_experiment_title()
        run_number = self.get_experiment_run_number()
        return f"{title} ({run_number})"

    def get_goniometer(self, idx: int = None) -> Goniometer:
        rotation_axis = (0.0, 1.0, 0.0)
        fixed_rotation = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        goniometer = GoniometerFactory.make_goniometer(rotation_axis, fixed_rotation)
        try:
            experiment_title = self.get_experiment_title()
            if "w=" in experiment_title:
                angle = float(experiment_title().split("w=")[1].split()[0])
            elif "wccr" in experiment_title:
                angle = float(experiment_title.split("wccr=")[1].split()[0])
            elif "wtl" in experiment_title:
                angle = float(experiment_title.split("wtl=")[1].split()[0])
            else:
                return goniometer
            goniometer.rotate_around_origin(rotation_axis, angle)
        except (ValueError, IndexError):
            pass
        return goniometer

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
            r, t = panel_projections[i + 1]
            r = tuple(map(int, r))
            t = tuple(map(int, t))
            panel.set_projection_2d(r, t)

        return detector

    def _get_pixel_size(self) -> Tuple[float, float]:
        # (mm)
        return (3.0, 3.0)

    def _get_image_size(self) -> Tuple[int, int]:
        # (px)
        return (64, 64)

    def _get_num_panels(self) -> int:
        return 11

    def _get_panel_names(self) -> List[str]:
        return ["%02d" % (i + 1) for i in range(11)]

    def _get_panel_gain(self):
        return 1.0

    def _get_panel_trusted_range(self) -> Tuple[int, int]:
        return (-1, 100000)

    def _get_panel_type(self) -> str:
        return "SENSOR_PAD"

    def _panel_0_params(self):
        if self._panel_0_flipped():
            return {
                "slow_axis": (0.793, 0.0, 0.609),
                "fast_axis": (0.0, -1.0, 0.0),
                "origin": (60.81, 96.0, -236.946),
            }
        return {
            "fast_axis": (-0.793, 0.0, -0.609),
            "slow_axis": (0.0, 1.0, 0.0),
            "origin": (213.099, -96.0, -120.041),
        }

    def _get_panel_origins(self) -> Tuple[Tuple[float, float, float]]:
        # (mm)
        return (
            self._panel_0_params()["origin"],
            (224.999, -96.0, 96.0),
            (60.809, -96.0, 236.945),
            (-214.172, -96.0, 118.198),
            (-224.999, -96.0, -96.0),
            (-60.809, -96.0, -236.945),
            (127.534, -256.614, 96.0),
            (-96.0, -256.614, 127.534),
            (-123.036, -258.801, -96.0),
            (96.0, -256.614, -127.534),
            (96.0, -278.0, 96.0),
        )

    def _panel_0_flipped(self) -> bool:
        if self._nxs_file["raw_data_1"]["run_number"][0] > 30000:
            return True
        return False

    def _get_panel_slow_axes(self) -> Tuple[Tuple[float, float, float]]:
        return (
            self._panel_0_params()["slow_axis"],
            (0.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.695, 0.719, -0.0),
            (0.0, 0.719, 0.695),
            (-0.707, 0.707, -0.0),
            (0.0, 0.719, -0.695),
            (-0.0, 0.0, -1.0),
        )

    def _get_panel_fast_axes(self) -> Tuple[Tuple[float, float, float]]:
        return (
            self._panel_0_params()["fast_axis"],
            (-0.0, -0.0, -1.0),
            (0.793, -0.0, -0.609),
            (0.788, -0.0, 0.616),
            (-0.0, -0.0, 1.0),
            (-0.793, -0.0, 0.609),
            (0.0, -0.0, -1.0),
            (1.0, -0.0, -0.0),
            (-0.0, -0.0, 1.0),
            (-1.0, -0.0, -0.0),
            (-1.0, -0.0, -0.0),
        )

    def _get_panel_projections_2d(self) -> dict:
        """
        Returns a projection of the
        detector flattened around the bottom panel (11)
        """

        p_w, p_h = self._get_image_size()
        panel_pos = {
            11: ((1, 0, 0, 1), (0, 0)),
            10: ((1, 0, 0, 1), (-p_h, 0)),
            8: ((1, 0, 0, 1), (p_h, 0)),
            7: ((1, 0, 0, 1), (0, -p_w)),
            9: ((1, 0, 0, 1), (0, p_w)),
            2: ((1, 0, 0, 1), (0, 2 * -p_w)),
            5: ((1, 0, 0, 1), (0, 2 * p_w)),
            3: ((1, 0, 0, 1), (p_h, 2 * -p_w)),
            4: ((1, 0, 0, 1), (p_h, 2 * p_w)),
            1: ((1, 0, 0, 1), (-p_h, 2 * -p_w)),
            6: ((1, 0, 0, 1), (-p_h, 2 * p_w)),
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

    def _get_sample_to_source_distance(self) -> float:
        # (mm)
        return 8300.0

    def _get_sample_to_source_direction(self) -> Tuple[float, float, float]:
        return (0, 0, -1)

    def _get_wavelength_range(self) -> Tuple[float, float]:
        # (A)
        return (0.2, 10.0)

    def get_scan(self, index=None) -> Scan:
        image_range = (1, self.get_num_images())
        properties = {"time_of_flight": self._get_time_of_flight()}
        return ScanFactory.make_scan_from_properties(
            image_range=image_range, properties=properties
        )

    def _get_time_channel_bins(self) -> List[float]:
        # (usec)
        return self._nxs_file["raw_data_1"]["instrument"]["dae"]["time_channels_1"][
            "time_of_flight"
        ][:]

    def _get_time_of_flight(self) -> Tuple[float]:
        # (usec)
        bins = self._get_time_channel_bins()
        return tuple(
            [float((bins[i] + bins[i + 1]) * 0.5) for i in range(len(bins) - 1)]
        )

    def get_num_images(self) -> int:
        return len(self._get_time_of_flight())

    def _load_raw_data(self):
        raw_data = []
        panel_size = self._get_image_size()
        total_pixels = panel_size[0] * panel_size[1]
        num_tof_bins = len(self._get_time_of_flight())

        # Panel positions are offset by 4 in raw_data array
        # See p24 of https://www.isis.stfc.ac.uk/Pages/sxd-user-guide6683.pdf
        idx_offset = 4

        for panel_idx in range(self._get_num_panels()):
            start_idx = (panel_idx * total_pixels) + (panel_idx * idx_offset)
            end_idx = start_idx + total_pixels
            panel_data = flex.int(
                self._nxs_file["raw_data_1/detector_1/counts"][0, start_idx:end_idx, :]
            )
            panel_data.reshape(flex.grid(panel_size[0], panel_size[1], num_tof_bins))
            raw_data.append(panel_data)
        self._raw_data = tuple(raw_data)

    def get_raw_data(self, index: int, use_loaded_data=True) -> Tuple[flex.int]:

        raw_data = []

        if use_loaded_data:
            if self._raw_data is None:
                self._load_raw_data()
            for panel in self._raw_data:
                data = panel[:, :, index : index + 1]
                data.reshape(flex.grid(panel.all()[0], panel.all()[1]))
                data.matrix_transpose_in_place()
                raw_data.append(data)

        else:
            panel_size = self._get_image_size()
            total_pixels = panel_size[0] * panel_size[1]

            # Panel positions are offset by 4 in raw_data array
            # See p24 of https://www.isis.stfc.ac.uk/Pages/sxd-user-guide6683.pdf
            idx_offset = 4

            for panel_idx in range(self._get_num_panels()):
                start_idx = (panel_idx * total_pixels) + (panel_idx * idx_offset)
                end_idx = start_idx + total_pixels
                panel_data = flex.int(
                    self._nxs_file["raw_data_1/detector_1/counts"][
                        0, start_idx:end_idx, index
                    ]
                )
                panel_data.reshape(flex.grid(panel_size[0], panel_size[1]))
                panel_data.matrix_transpose_in_place()
                raw_data.append(panel_data)

        return tuple(raw_data)
