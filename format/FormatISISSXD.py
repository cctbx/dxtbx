from __future__ import absolute_import, division, print_function

from sys import argv

import h5py
import numpy as np

from dials.array_family import flex

from dxtbx import IncorrectFormatError
from dxtbx.format.FormatNXTOFRAW import FormatNXTOFRAW
from dxtbx.model import Detector, Goniometer, Scan
from dxtbx.model.beam import BeamFactory


class FormatISISSXD(FormatNXTOFRAW):

    """
    Class to read NXTOFRAW files from the ISIS SXD
    (https://www.isis.stfc.ac.uk/Pages/sxd.aspx)

    """

    def __init__(self, image_file, **kwargs):
        super().__init__(image_file)
        if not FormatISISSXD.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self.nxs_file = self.open_file(image_file)
        self.detector = None
        self.raw_data = None

    def open_file(self, image_file):
        return h5py.File(image_file, "r")

    @staticmethod
    def understand(image_file):
        try:
            return FormatISISSXD.is_isissxd_file(image_file)
        except IOError:
            return False

    @staticmethod
    def is_isissxd_file(image_file):

        """
        Confirms if image_file is a NXTOFRAW format
        and from the SXD by confirming required fields
        are present and then checking the name attribute

        """

        def get_name(image_file):
            with h5py.File(image_file, "r") as handle:
                return handle["/raw_data_1/name"][0].decode()

        if not FormatNXTOFRAW.understand(image_file):
            return False

        return get_name(image_file) == "SXD"

    def load_raw_data(self):
        def get_detector_idx_array(detector_number, image_size, idx_offset):
            total_pixels = image_size[0] * image_size[1]
            min_range = (total_pixels * (detector_number - 1)) + (
                idx_offset * (detector_number - 1)
            )
            max_range = min_range + total_pixels
            return np.arange(min_range, max_range).reshape(image_size).T

        raw_counts = self.nxs_file["raw_data_1"]["detector_1"]["counts"][0, :, :]
        num_panels = self._get_num_panels()
        image_size = self._get_panel_size_in_px()

        # Index offset in SXD data
        # See p24 of https://www.isis.stfc.ac.uk/Pages/sxd-user-guide6683.pdf
        idx_offset = 4

        num_images = self.get_num_images()
        raw_data = []

        for i in range(1, num_panels + 1):
            idx_array = get_detector_idx_array(i, image_size, idx_offset)
            panel_array = np.zeros((idx_array.shape[0], idx_array.shape[1], num_images))
            for c_i, i in enumerate(idx_array):
                for c_j, j in enumerate(i):
                    panel_array[c_i, c_j, :] = raw_counts[j, :]
            flex_array = flex.double(np.ascontiguousarray(panel_array))
            flex_array.reshape(flex.grid(panel_array.shape))
            raw_data.append(flex_array)

        return tuple(raw_data)

    def get_raw_data(self, index):
        if self.raw_data is None:
            self.raw_data = self.load_raw_data()

        raw_data_idx = []
        for i in self.raw_data:
            arr = i[:, :, index : index + 1]
            arr.reshape(flex.grid(i.all()[0], i.all()[1]))
            raw_data_idx.append(arr)

        return tuple(raw_data_idx)

    def _get_detector(self):

        """
        Returns a  Detector instance with parameters taken from

        """

        num_panels = self._get_num_panels()
        panel_names = self._get_panel_names()
        panel_type = self._get_panel_type()
        image_size = self._get_panel_size_in_px()
        trusted_range = self._get_panel_trusted_range()
        pixel_size = self._get_panel_pixel_size_in_mm()
        fast_axes = self._get_panel_fast_axes()
        slow_axes = self._get_panel_slow_axes()
        panel_origins = self._get_panel_origins()
        gain = self._get_panel_gain()
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

        return detector

    """
    Hardcoded values not contained in the self.nxs_file are taken from
    https://doi.org/10.1107/S0021889806025921
    """

    def _get_time_channel_bins(self):
        return self.nxs_file["raw_data_1"]["instrument"]["dae"]["time_channels_1"][
            "time_of_flight"
        ][:]

    def _get_time_channels_in_seconds(self):
        bins = self._get_time_channel_bins()
        return [(bins[i] + bins[i + 1]) * 0.5 * 10 ** -6 for i in range(len(bins) - 1)]

    def _get_primary_flight_path_in_m(self):
        return 8.3

    def _get_num_panels(self):
        return 11

    def _get_panel_names(self):
        return ["%02d" % (i + 1) for i in range(11)]

    def _get_panel_origin_l2_vals_in_mm(self):
        return (
            262.787,
            262.787,
            262.787,
            262.787,
            262.787,
            302.212,
            302.212,
            302.212,
            302.212,
            302.212,
            311.178,
        )

    def _get_panel_gain(self):
        return 1.0

    def _get_panel_trusted_range(self):
        return (-1, 100000)

    def _get_panel_origins(self):
        return (
            (60.084, 103.128, -234.119),
            (222.316, 103.128, -94.855),
            (210.590, 103.128, 118.631),
            (-60.084, 103.128, 234.119),
            (-222.316, 103.128, 94.855),
            (-239.008, 101.244, -154.780),
            (257.166, -129.757, -91.437),
            (91.437, -129.757, 257.166),
            (-257.166, -129.757, 91.437),
            (-91.437, -129.757, -257.166),
            (32.732, -294.358, 95.467),
        )

    def _get_panel_slow_axes(self):
        return (
            (0.000, -1.000, 0.000),
            (0.000, -1.000, 0.000),
            (0.000, -1.000, 0.000),
            (0.000, -1.000, 0.000),
            (0.000, -1.000, 0.000),
            (0.000, -1.000, 0.000),
            (-0.707, -0.707, 0.000),
            (0.000, -0.707, -0.707),
            (0.707, -0.707, 0.000),
            (0.000, -0.707, 0.707),
            (0.000, 0.000, -1.000),
        )

    def _get_panel_fast_axes(self):
        return (
            (0.793, 0.000, 0.609),
            (0.000, 0.000, 1.000),
            (-0.793, 0.000, 0.609),
            (-0.793, 0.000, -0.609),
            (0.000, 0.000, -1.000),
            (0.793, 0.000, -0.609),
            (0.000, 0.000, 1.000),
            (-1.000, 0.000, 0.000),
            (0.000, 0.000, -1.000),
            (1.000, 0.000, 0.000),
            (-1.000, 0.000, 0.000),
        )

    def _get_s0(self):
        return (0, 0, 1)

    def _get_unit_s0(self):
        return tuple(self._get_s0() / np.linalg.norm(self._get_s0()))

    def _get_beam_direction(self):
        return (0, 0, 1)

    def _get_beam_polarization_normal(self):
        return (0, 0, 0)

    def _get_beam_polarization_fraction(self):
        return 0.5

    def _get_beam_flux(self):
        return 0.0

    def _get_beam_transmission(self):
        return 1.0

    def _get_beam_divergence(self):
        return 0.0

    def _get_beam_sigma_divergence(self):
        return 0.0

    def get_num_images(self):
        return len(self._get_time_channels_in_seconds())

    def get_beam(self, idx=None):
        s0 = self._get_s0()
        unit_s0 = self._get_unit_s0()
        wavelength = -1
        direction = self._get_beam_direction()
        divergence = self._get_beam_divergence()
        sigma_divergence = self._get_beam_sigma_divergence()
        polarization_normal = self._get_beam_polarization_normal()
        polarization_fraction = self._get_beam_polarization_fraction()
        flux = self._get_beam_flux()
        transmission = self._get_beam_transmission()

        beam = BeamFactory.make_beam(s0=s0, unit_s0=unit_s0, wavelength=wavelength)
        beam.set_direction(direction)
        beam.set_divergence(divergence)
        beam.set_sigma_divergence(sigma_divergence)
        beam.set_polarization_normal(polarization_normal)
        beam.set_polarization_fraction(polarization_fraction)
        beam.set_flux(flux)
        beam.set_transmission(transmission)

        return beam

    def get_detector(self, idx=None):
        return self._get_detector()

    def get_scan(self, idx=None):
        tof = self._get_time_channels_in_seconds()
        image_range = (1, len(tof))
        return Scan(tuple(map(int, image_range)), flex.double(list(map(float, tof))))

    def get_goniometer(self, idx=None):
        return Goniometer()

    def _get_panel_size_in_px(self):
        return (64, 64)

    def _get_panel_pixel_size_in_mm(self):
        return (3, 3)

    def _get_panel_type(self):
        return "SENSOR_PAD"


if __name__ == "__main__":
    for arg in argv[1:]:
        print(FormatISISSXD.understand(arg))
