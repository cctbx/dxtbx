from __future__ import annotations

import logging
from os.path import join

import h5py
import numpy as np

import cctbx.array_family.flex as flex

import dxtbx_flumpy as flumpy
from dxtbx import IncorrectFormatError
from dxtbx.format.FormatESSNMX import FormatESSNMX
from dxtbx.model.beam import BeamFactory, PolychromaticBeam, Probe
from dxtbx.model.scan import Scan, ScanFactory

logger = logging.getLogger(__name__)


class FormatESSNMXLaue(FormatESSNMX):
    """
    Class to read files from NMX
    https://europeanspallationsource.se/instruments/nmx
    flattened in the TOF dimension

    To add Laue data to a MANDI NeXus file:
    >>> from dxtbx.format.FormatMANDILaue import FormatMANDILaue
    >>> nxs_file = "/path/to/file/example.nxs.h5"
    >>> FormatESSNMXLaue.add_laue_data_to_nxs_file(nxs_file)

    """

    def __init__(self, image_file, **kwargs) -> None:
        if not FormatESSNMXLaue.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self._nxs_file = h5py.File(image_file, "r")
        self._raw_data = None

    @staticmethod
    def understand(image_file: str) -> bool:
        try:
            return FormatESSNMXLaue.is_nmx_laue_file(image_file)
        except (OSError, KeyError):
            return False

    @staticmethod
    def is_nmx_laue_file(image_file: str) -> bool:
        def get_name(image_file):
            try:
                with h5py.File(image_file, "r") as handle:
                    return handle["entry/instrument/name"][...].item().decode()
            except (OSError, KeyError, AttributeError):
                return ""

        if get_name(image_file) == "NMX":
            try:
                with h5py.File(image_file, "r") as handle:
                    return "laue_data" in handle["entry/instrument/detector_panel_0"]

            except (ValueError, IndexError):
                return False

    def get_num_images(self) -> int:
        return 1

    def get_scan(self, index=None) -> Scan:
        image_range = (1, self.get_num_images())
        return ScanFactory.make_scan_from_properties(
            image_range=image_range, properties={}
        )

    def get_beam(self, index: int | None = None) -> PolychromaticBeam:
        direction = self._get_sample_to_source_direction()
        distance = self._get_sample_to_source_distance()
        wavelength_range = self._get_wavelength_range()
        return BeamFactory.make_polychromatic_beam(
            direction=direction,
            sample_to_source_distance=distance,
            probe=Probe.neutron,
            wavelength_range=wavelength_range,
        )

    def get_raw_data(
        self,
        index: int,
    ) -> tuple[flex.int]:
        raw_data = []
        for panel in self._get_panels():
            spectra = panel["laue_data"][:]
            spectra = panel["laue_data"].astype(np.int32)
            raw_data.append(flumpy.from_numpy(np.ascontiguousarray(spectra)))

        return tuple(raw_data)

    @staticmethod
    def get_panel_names(nxs_file: h5py.File) -> list[str]:
        panel_names = []
        inst_dset = nxs_file["/entry/instrument/"]
        for name, dset in inst_dset.items():
            if dset.attrs.get("NX_class") == "NXdetector":
                panel_names.append(name)
        return panel_names

    @staticmethod
    def add_laue_data_to_nxs_file(
        nxs_file_path: str,
        output_name: str = "laue_data",
        compress=True,
    ) -> None:
        """
        Extracts event data and writes out laue data for each panel in
        nxs_file_path
        """

        print(f"Writing laue data in {nxs_file_path}")

        nxs_file = h5py.File(nxs_file_path, "r+")
        base_dir = "entry/instrument"

        panel_names = FormatESSNMXLaue.get_panel_names(nxs_file)
        for panel_name in panel_names:
            output_path = join(base_dir, panel_name)
            output_path = join(output_path, output_name)
            panel_data = FormatESSNMXLaue.generate_laue_data_for_panel(
                nxs_file, panel_name
            )
            if compress:
                nxs_file.create_dataset(
                    output_path, data=panel_data, compression="gzip"
                )
            else:
                nxs_file.create_dataset(output_path, data=panel_data)

        nxs_file.close()

    @staticmethod
    def generate_laue_data_for_panel(
        nxs_file: h5py.File,
        panel_name: str,
    ) -> np.array:
        """
        Generates laue data for a given panel
        """

        data = nxs_file[f"entry/instrument/{panel_name}/data"][:]

        laue_data = np.sum(data, axis=2)

        return laue_data
