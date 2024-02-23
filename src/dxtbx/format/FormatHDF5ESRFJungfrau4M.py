from __future__ import annotations

import sys

import h5py

from scitbx.array_family import flex

from dxtbx import flumpy
from dxtbx.format.FormatHDF5 import FormatHDF5


class FormatHDF5ESRFJungfrau4M(FormatHDF5):

    # A class to understand still-shot images from ESRF collected on a Jungfrau 4M.

    _cached_mask = None

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file, "r") as h5_handle:
            if len(h5_handle) != 1:
                return False
            key = list(h5_handle.keys())[0]
            if "instrument" not in h5_handle[key]:
                return False
            # instrument name is of form jungfrau4m_rr[X]_smx where X is empty, 4 or another number
            if not h5_handle[key]["instrument"].startswith("jungfrau4m_rr"):
                return False
            if not h5_handle[key]["instrument"].endswith("smx"):
                return False
        return True

    def _start(self):
        super()._start()
        image_file = self.get_image_file()
        self._h5_handle = h5py.File(image_file, "r")
        self.key = list(self._h5_handle.keys())[0]
        self.n_images = self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
            "data"
        ].shape[0]
        self.adus_per_photon = self._h5_handle[self.key]["instrument"][
            "jungfrau4m_rr_smx"
        ]["detector_information"]["adus_per_photon"]
        self.image_size = tuple(
            self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"]["data"].shape[
                1:
            ]
        )
        wavelength = self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
            "beam"
        ]["incident_wavelength"][()]
        x_pixel_size = (
            self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
                "detector_information"
            ]["x_pixel_size"][()]
            * 1000
        )  # convert m to mm
        y_pixel_size = (
            self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
                "detector_information"
            ]["y_pixel_size"][()]
            * 1000
        )  # convert m to mm
        distance = (
            self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
                "detector_information"
            ]["detector_distance"][()]
            * 1000
        )  # convert m to mm
        beam_center_x = self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
            "detector_information"
        ]["beam_center_x"][
            ()
        ]  # in px
        beam_center_y = self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
            "detector_information"
        ]["beam_center_y"][
            ()
        ]  # in px

        beam_center_x *= x_pixel_size
        beam_center_y *= y_pixel_size
        trusted_range = (
            self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
                "detector_information"
            ]["underload_value"][()],
            self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
                "detector_information"
            ]["saturation_value"][()],
        )
        exposure_time = self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
            "acquisition"
        ]["exposure_time"][()]

        self._detector_model = self._detector_factory.simple(
            sensor="UNKNOWN",
            distance=distance,
            beam_centre=(
                beam_center_x,
                beam_center_y,
            ),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(
                x_pixel_size,
                y_pixel_size,
            ),
            image_size=(self.image_size[1], self.image_size[0]),
            trusted_range=trusted_range,
            mask=self.get_static_mask(),
        )
        self._beam_model = self._beam_factory.simple(wavelength)
        self._scan_model = self._scan_factory.make_scan(
            image_range=(1, self.n_images),
            exposure_times=exposure_time,
            oscillation=(0.0, 0.0),
            epochs=list(range(self.n_images)),
        )
        self._goniometer_model = self._goniometer_factory.known_axis((0, 1, 0))

    def get_raw_data(self, index=None):
        if index is None:
            index = 0
        # data can be int32 with adus_per_photon != 1.0 or float16 with adus_per_photon == 1.0
        data = (
            (
                self._h5_handle[self.key]["measurement"]["data"][index]
                / self.adus_per_photon
            )
            if self.adus_per_photon != 1.0
            else self._h5_handle[self.key]["measurement"]["data"][index]
        )
        return flex.double(data.astype(float))

    def get_num_images(self):
        return self.n_images

    def get_beam(self, index=None):
        return self._beam(index)

    def _beam(self, index=None):
        return self._beam_model

    def get_detector(self, index=None):
        return self._detector(index)

    def get_static_mask(self):
        if FormatHDF5ESRFJungfrau4M._cached_mask is None:
            mask = self._h5_handle[self.key]["instrument"]["jungfrau4m_rr_smx"][
                "detector_information"
            ]["pixel_mask"]
            mask = flumpy.from_numpy(mask[()])
            mask_array = mask == 0
            FormatHDF5ESRFJungfrau4M._cached_mask = mask_array
        return FormatHDF5ESRFJungfrau4M._cached_mask

    def _detector(self, index=None):
        return self._detector_model

    def _goniometer(self):
        return self._goniometer_model

    def _scan(self):
        return self._scan_model


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatHDF5ESRFJungfrau4M.understand(arg))
