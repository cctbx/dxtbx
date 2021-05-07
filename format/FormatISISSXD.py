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

    def load_raw_data(self, as_numpy_arrays=False):
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

        for n in range(1, num_panels + 1):
            idx_array = get_detector_idx_array(n, image_size, idx_offset)
            panel_array = np.zeros((idx_array.shape[0], idx_array.shape[1], num_images))
            for c_i, i in enumerate(idx_array):
                for c_j, j in enumerate(i):
                    panel_array[c_i, c_j, :] = raw_counts[j, :]
            if as_numpy_arrays:
                raw_data.append(panel_array)
            else:
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
            (60.843, -96.0, -236.969),
            (225.0, 96.0, 96.0),
            (60.843, 96.0, 236.969),
            (-214.172, 96.0, 118.166),
            (-226.0, 96.0, -96.0),
            (-60.843, 96.0, -236.969),
            (127.502, 256.582, 96.0),
            (-93.0, 256.582, 127.502),
            (-124.047, 258.791, -96.0),
            (94.0, 256.582, -127.502),
            (95.0, 279.0, 96.0),
        )

    def _get_panel_slow_axes(self):
        return (
            (0.793, 0.0, 0.609),
            (0.0, -1.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.695, -0.719, -0.0),
            (0.0, -0.719, 0.695),
            (-0.707, -0.707, -0.0),
            (0.0, -0.719, -0.695),
            (-0.0, 0.0, -1.0),
        )

    def _get_panel_fast_axes(self):
        return (
            (0.0, 1.0, 0.0),
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

    def _get_s0(self):
        return (0, 0, 0)

    def _get_unit_s0(self):
        return (0, 0, 1)
        return tuple(self._get_s0() / np.linalg.norm(self._get_s0()))

    def _get_beam_direction(self):
        return (0, 0, -1)

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
        wavelength = 1
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

    def _get_raw_spectra_array(self):
        # Returns 2D array of (pixels, time_channels) for all 11 detectors
        return self.nxs_file["raw_data_1"]["detector_1"]["counts"][:][0]

    def _get_panel_images(self):

        """
        Returns a list of arrays (x_num_px, y_num_px, num_time_channels)
        for each panel, ordered from 1-11
        """
        raw_data = self._get_raw_spectra_array()

        # Panel positions are offset by 4 in raw_data array
        # See p24 of https://www.isis.stfc.ac.uk/Pages/sxd-user-guide6683.pdf
        panel_size = self._get_panel_size_in_px()
        total_px = panel_size[0] * panel_size[1]
        offsets = [
            (((total_px * i) + (i * 4)), ((total_px * (i + 1)) + (i * 4)))
            for i in range(11)
        ]
        panel_raw_data = [raw_data[i[0] : i[1], :] for i in offsets]

        panel_size = self._get_panel_size_in_px()
        time_channel_size = len(self._get_time_channels_in_seconds())
        array_shape = (panel_size[0], panel_size[1], time_channel_size)
        return [i.reshape(array_shape) for i in panel_raw_data]

    def get_reflection_table_from_use_file(self, use_file, specific_panel=None):

        import dials_array_family_flex_ext
        from scipy import interpolate

        import cctbx.array_family.flex

        def is_data_row(row):
            num_data_columns = 13
            return len(row.split()) == num_data_columns

        def is_active_row(count, lines):
            num_active_row_columns = 9
            active_pos = -1
            if len(lines[count + 1].split()) == num_active_row_columns:
                return lines[count + 1].split()[active_pos] == "1"
            return False

        def get_hkl(row):
            hkl_pos = (1, 4)
            assert is_data_row(row), "Cannot extract data from this row"
            hkl = row.split()[hkl_pos[0] : hkl_pos[1]]
            return tuple(map(int, hkl))

        def get_single_float_value(row, idx):
            # assert(is_data_row(row)), "Cannot extract data from this row"
            return float(row.split()[idx])

        def get_x(row):
            x_pos = 4
            return get_single_float_value(row, x_pos)

        def get_dx(row):
            dx_pos = 3
            return get_single_float_value(row, dx_pos)

        def get_y(row):
            y_pos = 5
            return get_single_float_value(row, y_pos)

        def get_dy(row):
            dy_pos = 4
            return get_single_float_value(row, dy_pos)

        def get_tof_curve_coefficients(tof_vals):
            x = [i + 1 for i in range(len(tof_vals))]
            return interpolate.splrep(tof_vals, x)

        def get_tof_frame(tof, tof_curve_coeffs):
            return float(interpolate.splev(tof, tof_curve_coeffs))

        def get_bbox(x, dx, y, dy, tof, dtof, tof_curve_coeffs):
            frame = get_tof_frame(tof, tof_curve_coeffs)
            dframe1 = get_tof_frame(tof - dtof, tof_curve_coeffs)
            dframe2 = get_tof_frame(tof + dtof, tof_curve_coeffs)
            dframe = dframe2 - dframe1
            return (
                int(x - dx),
                int(x + dx),
                int(y - dy),
                int(y + dy),
                int(frame - dframe),
                int(frame + dframe),
            )

        def get_tof(row):
            tof_pos = 6
            return get_single_float_value(row, tof_pos)

        def get_dtof(row):
            dtof_pos = 5
            return get_single_float_value(row, dtof_pos)

        def get_tof_wavelength_in_ang(L0, L, tof):
            h = 6.62607004e-34
            m_n = 1.67493e-27
            return ((h * tof) / (m_n * (L0 + L))) * 10 ** 10

        def get_pixel_wavelength_in_ang(
            x, y, tof, L0, centroid_l, pixel_size_in_mm, panel_size_in_px
        ):
            import numpy as np

            rel_x = abs(x) * pixel_size_in_mm[0] * 10 ** -3
            rel_y = abs(y) * pixel_size_in_mm[1] * 10 ** -3
            rel_pos = np.sqrt(np.square(rel_x) + np.square(rel_y))
            rel_L = np.sqrt(np.square(rel_pos) + np.square(centroid_l))

            return get_tof_wavelength_in_ang(L0, rel_L, tof)

        def get_panel_l(panel_idx):
            panel_l_vals = [
                0.225,
                0.225,
                0.225,
                0.225,
                0.225,
                0.225,
                0.270,
                0.270,
                0.270,
                0.270,
                0.278,
            ]
            return panel_l_vals[panel_idx]

        def at_start_of_table(row):
            return row.startswith(" NSEQ  ")

        def at_end_of_table(row):
            return row.startswith("DETECTOR CALIBRATION FILE ")

        def convert_coord_to_dials(x, y):
            offset = 32
            return (x + offset, y + offset)

        def get_data(use_file, specific_panel):

            with open(use_file, "r") as g:
                lines = g.readlines()

            panel_size = (64, 64)
            data_names = [
                "x",
                "y",
                "frame",
                "tof",
                "panel",
                "hkl",
                "wavelength",
                "bbox",
            ]
            data = {i: [] for i in data_names}

            recording_table = False
            panel_num = -1
            tof_vals = self._get_time_channels_in_seconds()
            tof_vals = [i * 10 ** 6 for i in tof_vals]
            tof_curve_coeffs = get_tof_curve_coefficients(tof_vals)

            for count, line in enumerate(lines):

                if recording_table:
                    if is_data_row(line) and is_active_row(count, lines):
                        if specific_panel is not None:
                            if specific_panel != panel_num:
                                continue
                        data["hkl"].append(get_hkl(line))
                        if panel_num == 0:
                            y = -get_x(line)
                            x = -get_y(line)
                            dy = get_dx(lines[count + 1])
                            dx = get_dy(lines[count + 1])
                        else:
                            x = get_x(line)
                            y = get_y(line)
                            dx = get_dx(lines[count + 1])
                            dy = get_dy(lines[count + 1])

                        tof = get_tof(line)
                        wavelength = get_pixel_wavelength_in_ang(
                            x,
                            y,
                            tof * 10 ** -6,
                            8.3,
                            get_panel_l(panel_num),
                            pixel_size,
                            panel_size,
                        )
                        dtof = get_dtof(lines[count + 1])
                        x, y = convert_coord_to_dials(x, y)
                        bbox = get_bbox(x, dx, y, dy, tof, dtof, tof_curve_coeffs)

                        data["panel"].append(panel_num)
                        data["tof"].append(tof)
                        data["wavelength"].append(wavelength)
                        data["x"].append(x)
                        data["y"].append(y)
                        data["frame"].append(get_tof_frame(tof, tof_curve_coeffs))
                        data["bbox"].append(bbox)

                if at_start_of_table(line):
                    recording_table = True

                if at_end_of_table(line):
                    panel_num += 1
                    recording_table = False

            return data

        def calc_s0(unit_s0, wavelength):
            return unit_s0 * (1.0 / wavelength)

        pixel_size = (3, 3)
        data = get_data(use_file, specific_panel)
        nrows = len(data["hkl"])

        unit_s0 = np.array(self._get_unit_s0())
        use_reflection_table = (
            dials_array_family_flex_ext.reflection_table.empty_standard(nrows)
        )
        use_reflection_table["xyzobs.px.value"] = cctbx.array_family.flex.vec3_double(
            nrows
        )
        use_reflection_table["xyzobs.mm.value"] = cctbx.array_family.flex.vec3_double(
            nrows
        )
        use_reflection_table["tof_wavelength"] = cctbx.array_family.flex.double(nrows)
        use_reflection_table["tof"] = cctbx.array_family.flex.double(nrows)
        use_reflection_table["tof_s0"] = cctbx.array_family.flex.vec3_double(nrows)
        use_reflection_table["bbox"] = dials_array_family_flex_ext.int6(nrows)
        use_reflection_table["miller_index"] = cctbx.array_family.flex.miller_index(
            nrows
        )
        use_reflection_table["panel"] = cctbx.array_family.flex.size_t(nrows)
        use_reflection_table["flags"] = cctbx.array_family.flex.size_t(nrows, 32)

        for i in range(nrows):
            use_reflection_table["xyzobs.px.value"][i] = (
                data["x"][i],
                data["y"][i],
                data["frame"][i],
            )
            use_reflection_table["xyzobs.mm.value"][i] = (
                data["x"][i] * pixel_size[0],
                data["y"][i] * pixel_size[1],
                data["tof"][i],
            )
            use_reflection_table["tof_wavelength"][i] = data["wavelength"][i]
            use_reflection_table["tof"][i] = data["tof"][i]
            use_reflection_table["tof_s0"][i] = calc_s0(unit_s0, data["wavelength"][i])
            use_reflection_table["miller_index"][i] = data["hkl"][i]
            use_reflection_table["panel"][i] = data["panel"][i]
            use_reflection_table["bbox"][i] = data["bbox"][i]

        use_reflection_table.set_flags(
            use_reflection_table["miller_index"] != (0, 0, 0),
            use_reflection_table.flags.indexed,
        )

        return use_reflection_table


if __name__ == "__main__":
    for arg in argv[1:]:
        print(FormatISISSXD.understand(arg))
