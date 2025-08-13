from __future__ import annotations

import xml
import xml.etree.ElementTree as ET
from multiprocessing import Pool, cpu_count
from os.path import join
from sys import argv

import h5py
import numpy as np

from cctbx.array_family import flex

import dxtbx_flumpy as flumpy
from dxtbx import IncorrectFormatError
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.model import Detector
from dxtbx.model.beam import BeamFactory, PolychromaticBeam, Probe
from dxtbx.model.goniometer import Goniometer, GoniometerFactory
from dxtbx.model.scan import Scan, ScanFactory
from dxtbx.model.tof_helpers import InstrumentDefinitionReader


class FormatMANDI(FormatHDF5):
    """
    Class to read NXTOFRAW from MaNDi
    (https://neutrons.ornl.gov/mandi)
    """

    def __init__(self, image_file: str, **kwargs) -> None:
        super().__init__(image_file)
        if not FormatMANDI.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self.image_file = image_file
        self.nxs_file = self.open_file(image_file)
        self._base_entry = self.get_base_entry_name(self.nxs_file)
        self.detector = None
        self.raw_data = None
        self.xml_reader = InstrumentDefinitionReader()
        self.xml_file = self.get_xml_file()

    def open_file(self, image_file_path: str) -> h5py.File:
        return h5py.File(image_file_path, "r")

    def get_base_entry_name(self, nxs_file: h5py.File) -> str:
        return list(nxs_file.keys())[0]

    def get_xml_file(self) -> xml.etree.ElementTree.Element:
        xml_string = self.nxs_file[self._base_entry]["instrument/instrument_xml/data"][
            0
        ].decode()
        return ET.fromstring(xml_string)

    @staticmethod
    def understand(image_file: str) -> bool:
        try:
            return FormatMANDI.is_mandi_file(image_file)
        except OSError:
            return False

    @staticmethod
    def is_mandi_file(image_file: str) -> bool:
        def get_name(image_file: str) -> str:
            with h5py.File(image_file, "r") as handle:
                if len(handle) == 0:
                    return ""
                base_entry = list(handle.keys())[0]
                if f"{base_entry}/instrument/name" not in handle:
                    return ""
                try:
                    return handle[f"{base_entry}/instrument/name"][0].decode()
                except (ValueError, IndexError):
                    return ""

        return get_name(image_file) == "MANDI"

    def get_instrument_name(self) -> str:
        return self.nxs_file[self._base_entry]["instrument/name"][0].decode()

    def get_experiment_title(self) -> str:
        return self.nxs_file[self._base_entry]["title"][0].decode()

    def get_experiment_run_number(self) -> str:
        return str(self.nxs_file[self._base_entry]["run_number"][0].decode())

    def get_experiment_description(self) -> str:
        title = self.get_experiment_title()
        run_number = self.get_experiment_run_number()
        return f"{title} ({run_number})"

    def get_raw_data(self, index: int) -> tuple[flex.int]:
        raw_data = []
        panel_size = self._get_image_size()
        for panel_name in self._get_panel_names():
            spectra = np.reshape(
                self.nxs_file[self._base_entry][f"{panel_name}_events"]["spectra"][
                    :, index : index + 1
                ],
                panel_size,
            ).T
            raw_data.append(flumpy.from_numpy(np.ascontiguousarray(spectra)))

        return tuple(raw_data)

    def get_detector(self, index: int = None) -> Detector:
        num_panels = self._get_num_panels()
        panel_names = self._get_panel_names()
        panel_type = self._get_panel_type()
        image_size = self._get_image_size()
        trusted_range = self._get_panel_trusted_range()
        pixel_size = self._get_pixel_size()
        panel_origins, fast_axes, slow_axes = (
            self.xml_reader.get_dials_detector_geometry(self.xml_file)
        )
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

    def _get_time_channel_bins(self) -> list[float]:
        # (usec)
        return self.nxs_file[self._base_entry]["time_of_flight"][:]

    def _get_time_of_flight(self) -> list[float]:
        # (usec)
        bins = self._get_time_channel_bins()
        return tuple([float(bins[i] + bins[i + 1]) * 0.5 for i in range(len(bins) - 1)])

    def _get_sample_to_source_distance(self) -> float:
        # (mm)
        return 30000

    def _get_sample_to_source_direction(self) -> tuple[float, float, float]:
        return (0, 0, -1)

    def _get_num_panels(self) -> int:
        return self.xml_reader.get_num_panels(self.xml_file)

    def _get_panel_names(self) -> tuple[str]:
        return self.xml_reader.get_panel_names(self.xml_file)

    def _get_panel_gain(self) -> float:
        return 1.0

    def _get_panel_trusted_range(self) -> tuple[int, int]:
        return (-1, 100000)

    def get_num_images(self) -> int:
        return len(self._get_time_of_flight())

    def get_beam(self, index: int = None) -> PolychromaticBeam:
        direction = self._get_sample_to_source_direction()
        distance = self._get_sample_to_source_distance()
        wavelength_range = (2.0, 4.0)
        return BeamFactory.make_polychromatic_beam(
            direction=direction,
            wavelength_range=wavelength_range,
            sample_to_source_distance=distance,
            probe=Probe.neutron,
        )

    def get_scan(self, index=None) -> Scan:
        image_range = (1, self.get_num_images())
        properties = {"time_of_flight": self._get_time_of_flight()}
        return ScanFactory.make_scan_from_properties(
            image_range=image_range, properties=properties
        )

    def get_goniometer(self, idx: int = None) -> Goniometer:
        rotation_axis_phi = (0.0, 1.0, 0.0)
        rotation_axis_omega = (0.0, 1.0, 0.0)
        rotation_axis_chi = (0.0, 0.0, 1.0)
        fixed_rotation = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

        goniometer = GoniometerFactory.make_goniometer(
            rotation_axis_phi, fixed_rotation
        )
        phi = self.nxs_file["entry/DASlogs/phi/average_value"][0]
        omega = self.nxs_file["entry/DASlogs/omega/average_value"][0]
        chi = self.nxs_file["entry/DASlogs/chi/average_value"][0]

        goniometer.rotate_around_origin(rotation_axis_phi, phi)
        goniometer.rotate_around_origin(rotation_axis_chi, chi)
        goniometer.rotate_around_origin(rotation_axis_omega, omega)
        return goniometer

    def _get_image_size(self) -> tuple[int, int]:
        return (256, 256)

    def _get_pixel_size(self) -> tuple[float, float]:
        return (0.618, 0.618)

    def _get_panel_type(self) -> str:
        return "SENSOR_PAD"

    def _get_panel_projections_2d(self) -> dict:
        p_w, p_h = self._get_image_size()
        panel_pos = {}
        count = 1
        for i in range(8):
            for j in range(5):
                panel_pos[count] = ((1, 0, 0, 1), (p_h * i, p_w * j))
                count += 1
        panel_pos[count] = ((1, 0, 0, 1), (p_h * i, p_w * j))
        return panel_pos

    @staticmethod
    def add_histogram_data_to_nxs_file(
        nxs_file_path: str,
        remove_event_data: bool,
        spectra_output_name: str = "spectra",
        write_tof_bins: bool = True,
        delta_tof: float = 5,  # (usec)
        tof_padding: float = 100,  # (usec)
        panel_size: tuple[int, int] = (256, 256),  # (px)
        nproc: int = 8,
    ) -> None:
        tof_bins = FormatMANDI.generate_tof_bins(
            nxs_file=nxs_file_path,
            panel_size=panel_size,
            delta_tof=delta_tof,
            padding=tof_padding,
        )
        FormatMANDI.write_histogram_data(
            nxs_file_path=nxs_file_path,
            tof_bins=tof_bins,
            panel_size=panel_size,
            remove_event_data=remove_event_data,
            spectra_output_name=spectra_output_name,
            write_tof_bins=write_tof_bins,
        )

    @staticmethod
    def write_histogram_data(
        nxs_file_path: str,
        tof_bins: np.array,
        panel_size: tuple[int, int],
        remove_event_data: bool,
        spectra_output_name: str = "spectra",
        write_tof_bins: bool = True,
        nproc: int = 8,
    ) -> None:
        """
        Generates histogram spectra for a given detector and writes it to nxs_file
        """

        def delete_event_data(nxs_file, base_dir, panel_name):
            del nxs_file[join(join(base_dir, panel_name), "event_index")]
            del nxs_file[join(join(base_dir, panel_name), "event_id")]
            del nxs_file[join(join(base_dir, panel_name), "event_time_zero")]
            del nxs_file[join(join(base_dir, panel_name), "event_time_offset")]

        print(f"Writing histrogram data in {nxs_file_path}")
        print(f"Remove event data: {remove_event_data}")
        nxs_file = h5py.File(nxs_file_path, "r+")
        base_dir = list(nxs_file.keys())[0]

        panel_names = FormatMANDI.get_panel_names(nxs_file)
        written_tof_bins = False
        for panel_name in panel_names:
            print(f"Processing panel {panel_name}")
            output_path = join(base_dir, panel_name)
            output_path = join(output_path, spectra_output_name)
            print(f"Writing spectra to {output_path}")
            panel_spectra = FormatMANDI.generate_histogram_data_for_panel(
                nxs_file, tof_bins, panel_size, panel_name, nproc
            )
            nxs_file.create_dataset(output_path, data=panel_spectra, compression="gzip")
            if remove_event_data:
                delete_event_data(nxs_file, base_dir, panel_name)
                print(f"Removed event data for {panel_name}")
            if write_tof_bins and not written_tof_bins:
                tof_path = join(base_dir, "time_of_flight")
                print(f"Writing time of flight bins to {tof_path}")
                nxs_file.create_dataset(tof_path, data=tof_bins, compression="gzip")
                written_tof_bins = True

        nxs_file.close()

    @staticmethod
    def compute_event_histogram(
        args: tuple[int, np.array, np.array, np.array],
    ) -> np.array:
        pixel_idx, event_time_offset, corrected_event_id, tof_bins = args
        h, _ = np.histogram(
            event_time_offset[corrected_event_id == pixel_idx], tof_bins
        )
        return h

    @staticmethod
    def generate_histogram_data_for_panel(
        nxs_file: h5py.File,
        tof_bins: np.array,
        panel_size: tuple[int, int],
        panel_name: str,
        nproc=8,
    ) -> np.array:
        """
        Generates histogram data for a given panel
        """

        ## Get panel data
        panel_number = FormatMANDI.get_panel_number(panel_name)
        # Actual pixel ids, starting from bottom left and going up y axis first
        event_id = nxs_file[f"entry/{panel_name}/event_id"][:]
        # Time each event_id was triggered after event_time_zero (sec)
        event_time_offset = nxs_file[f"entry/{panel_name}/event_time_offset"][:]

        # event_ids are given with an offset
        event_id_offset = panel_number * panel_size[0] * panel_size[1]
        corrected_event_id = event_id - event_id_offset

        num_pixels = panel_size[0] * panel_size[1]
        spectra = np.zeros((num_pixels, len(tof_bins) - 1), dtype=np.int32)

        num_cpu = cpu_count()
        if nproc > num_cpu:
            nproc = num_cpu

        pool = Pool(processes=nproc)

        args_list = [
            (i, event_time_offset, corrected_event_id, tof_bins)
            for i in range(num_pixels)
        ]
        results = pool.map(FormatMANDI.compute_event_histogram, args_list)

        pool.close()
        pool.join()

        for i, h in enumerate(results):
            spectra[i] = h

        return spectra

    @staticmethod
    def get_time_range_for_panel(
        nxs_file: h5py.File, panel_size: tuple[float, float], panel_name: str
    ) -> tuple[float, float]:
        """
        Returns the range of event times for a given panel
        """

        def event_data_is_valid(event_id, event_time_offset):
            if len(event_id) == 0 or len(event_time_offset) == 0:
                return False
            return len(event_id) == len(event_time_offset)

        panel_number = FormatMANDI.get_panel_number(panel_name)
        event_index = nxs_file[f"entry/{panel_name}/event_index"]
        event_id = nxs_file[f"entry/{panel_name}/event_id"]
        event_time_zero = nxs_file[f"entry/{panel_name}/event_time_zero"]
        event_time_offset = nxs_file[f"entry/{panel_name}/event_time_offset"]

        if not event_data_is_valid(event_id, event_time_offset):
            return None, None

        num_pixels = panel_size[0] * panel_size[1]
        event_id_offset = panel_number * num_pixels - 1

        raw_event_id = event_id[event_index[0]]
        corrected_event_id = raw_event_id - event_id_offset
        min_event_time = event_time_zero[0] + event_time_offset[corrected_event_id]

        max_idx = int(event_index[-1] - 1)
        raw_event_id = event_id[max_idx]
        corrected_event_id = raw_event_id - event_id_offset
        max_event_time = (
            event_time_zero[max_idx] + event_time_offset[corrected_event_id]
        )

        return min_event_time, max_event_time

    @staticmethod
    def get_time_range_for_dataset(
        nxs_file_path: str, panel_size: tuple[int, int]
    ) -> tuple[float, float]:
        """
        Iterates over num_panels to find the overall min/max tof event recorded
        """

        # TODO get panel_size and panel_number from nxs_file xml

        nxs_file = h5py.File(nxs_file_path, "r")

        min_tof = -1
        max_tof = -1

        panel_names = FormatMANDI.get_panel_names(nxs_file)

        for panel_name in panel_names:
            try:
                min_event_time, max_event_time = FormatMANDI.get_time_range_for_panel(
                    nxs_file, panel_size, panel_name
                )
                if min_event_time is None or max_event_time is None:
                    # Some banks contain no data
                    continue
                if min_tof < 0 or min_event_time < min_tof:
                    min_tof = min_event_time
                if max_event_time > max_tof:
                    max_tof = max_event_time
            except KeyError:
                # Some banks not always present
                pass

        nxs_file.close()

        return min_tof, max_tof

    @staticmethod
    def generate_tof_bins(
        nxs_file: str,
        panel_size: tuple[float, float],
        delta_tof: float = 50,
        padding: float = 100,
    ) -> np.ndarray:
        """
        delta_tof: float (usec)
        padding: float (usec)
        """

        min_tof, max_tof = FormatMANDI.get_time_range_for_dataset(nxs_file, panel_size)
        min_tof = min_tof - padding
        max_tof = max_tof + padding
        print(
            f"Time of flight range for {nxs_file}: {round(min_tof, 3)} - {round(max_tof, 3)} (usec)"
        )
        num_bins = int((max_tof - min_tof) / delta_tof)
        return np.linspace(min_tof, max_tof, num_bins)

    @staticmethod
    def get_panel_names(nxs_file: h5py.File) -> list[str]:
        raw_names = [i for i in nxs_file[list(nxs_file.keys())[0]] if "bank" in i]
        names = []
        for name in raw_names:
            try:
                entry = (name, int("".join([j for j in name if j.isdigit()])))
                names.append(entry)
            except ValueError:  # Other fields with bank in them
                continue
        return [i[0] for i in sorted(names, key=lambda x: x[0])]

    @staticmethod
    def get_panel_number(panel_name: str) -> int:
        return int("".join([i for i in panel_name if i.isdigit()]))


if __name__ == "__main__":
    for arg in argv[1:]:
        print(FormatMANDI.understand(arg))
