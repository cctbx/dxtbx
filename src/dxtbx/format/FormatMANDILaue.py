from __future__ import annotations

from os.path import join

import h5py
import numpy as np
from dials.array_family import flex

from dxtbx import flumpy
from dxtbx.format.FormatMANDI import FormatMANDI
from dxtbx.model.beam import BeamFactory, PolychromaticBeam, Probe
from dxtbx.model.scan import Scan, ScanFactory


class FormatMANDILaue(FormatMANDI):
    """
    Class to read MANDI data flattened along the ToF dimension
    to approximate data expected from IMAGINE-X
    https://neutrons.ornl.gov/imagine

    To add Laue data to a MANDI NeXus file:
    >>> from dxtbx.format.FormatMANDILaue import FormatMANDILaue
    >>> nxs_file = "/path/to/file/example.nxs.h5"
    >>> FormatMANDILaue.add_laue_data_to_nxs_file(nxs_file)
    """

    @staticmethod
    def understand(image_file: str) -> bool:
        try:
            return FormatMANDILaue.is_mandi_laue_file(image_file)
        except OSError:
            return False

    @staticmethod
    def is_mandi_laue_file(image_file: str) -> bool:
        with h5py.File(image_file, "r") as handle:
            # File is not empty
            if len(handle) == 0:
                return False
            base_entry = list(handle.keys())[0]
            # Has a name entry
            if f"{base_entry}/instrument/name" not in handle:
                return False
            # Name can be decoded as is MANDI
            try:
                name = handle[f"{base_entry}/instrument/name"][0].decode()
                if name == "MANDI" and "bank10_events" in handle[base_entry].keys():
                    # Laue data present
                    return "laue_data" in handle[f"{base_entry}/bank10_events"]

            except (ValueError, IndexError):
                return False

    def get_raw_data(self, index: int) -> tuple[flex.int]:
        raw_data = []
        for panel_name in self._get_panel_names():
            panel_data = self.nxs_file[self._base_entry][f"{panel_name}_events"][
                "laue_data"
            ][:]
            panel_data = panel_data.astype(np.int32)
            raw_data.append(flumpy.from_numpy(np.ascontiguousarray(panel_data)))
        return tuple(raw_data)

    def get_scan(self, index: int | None = None) -> Scan:
        image_range = (1, self.get_num_images())
        return ScanFactory.make_scan_from_properties(
            image_range=image_range, properties={}
        )

    def get_num_images(self) -> int:
        return 1

    def get_beam(self, index: int | None = None) -> PolychromaticBeam:
        direction = self._get_sample_to_source_direction()
        distance = self._get_sample_to_source_distance()
        wavelength_range = (2.0, 4.0)
        return BeamFactory.make_polychromatic_beam(
            direction=direction,
            sample_to_source_distance=distance,
            probe=Probe.neutron,
            wavelength_range=wavelength_range,
        )

    @staticmethod
    def add_laue_data_to_nxs_file(
        nxs_file_path: str,
        output_name: str = "laue_data",
        panel_size: tuple[int, int] = (256, 256),  # (px)
        compress=True,
    ) -> None:
        """
        Extracts event data and writes out laue data for each panel in
        nxs_file_path
        """

        print(f"Writing laue data in {nxs_file_path}")

        nxs_file = h5py.File(nxs_file_path, "r+")
        base_dir = list(nxs_file.keys())[0]

        panel_names = FormatMANDI.get_panel_names(nxs_file)
        for panel_name in panel_names:
            output_path = join(base_dir, panel_name)
            output_path = join(output_path, output_name)
            panel_data = FormatMANDI.generate_laue_data_for_panel(
                nxs_file, panel_size, panel_name
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
        panel_size: tuple[int, int],
        panel_name: str,
    ) -> np.array:
        """
        Generates laue data for a given panel
        """

        ## Get panel data
        panel_number = FormatMANDI.get_panel_number(panel_name)
        # Actual pixel ids, starting from bottom left and going up y axis first
        event_id = nxs_file[f"entry/{panel_name}/event_id"][:]

        # event_ids are given with an offset
        event_id_offset = panel_number * panel_size[0] * panel_size[1]
        corrected_event_id = event_id - event_id_offset

        num_pixels = panel_size[0] * panel_size[1]

        counts = np.bincount(corrected_event_id)

        # Pad if no events for pixels outside of max event id
        if len(counts) < num_pixels:
            counts = np.pad(counts, (0, num_pixels - len(counts)), mode="constant")
        else:
            counts = counts[:num_pixels]

        laue_array = np.reshape(counts, panel_size)

        return laue_array
