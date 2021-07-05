"""An implementation of the CBF image reader for Pilatus images, for the P12M-DLS"""


import binascii
import calendar
import math
import sys

import libtbx
from cctbx.eltbx import attenuation_coefficient
from scitbx import matrix

from dxtbx.ext import uncompress
from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp
from dxtbx.masking import GoniometerMaskerFactory
from dxtbx.model import Detector, ParallaxCorrectedPxMmStrategy


class FormatCBFMiniPilatusDLS12M(FormatCBFMiniPilatus):
    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 120-0100" in header
            ):
                return True

        return False

    @staticmethod
    def has_dynamic_shadowing(**kwargs):
        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if dynamic_shadowing in (libtbx.Auto, "Auto"):
            return True
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        # if multi_panel == False, then interpret data as 24 panels, where each
        # row of 5 panels is grouped as one "panel"
        # elif multi_panel == True, then interpret data as 120 panels,
        # 24 rows * 5 columns
        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)
        self._multi_panel = kwargs.get("multi_panel", False)

        super().__init__(image_file, **kwargs)

    def _detector(self):

        # module positions from detector blueprints - modelling at the moment as
        # 24 modules, each consisting of 5 sensors (the latter is ignored)
        x = matrix.col((-1, 0, 0))
        y = matrix.col((0, 1, 0))
        z = matrix.col((0, 0, 1))

        beam_xy = self._cif_header_dictionary["Beam_xy"]
        beam_xy = beam_xy.replace("(", "").replace(")", "").replace(",", "").split()[:2]
        obs_beam_x, obs_beam_y = [float(f) for f in beam_xy]

        ideal_beam_x = 1075
        ideal_beam_y = 2594

        beam_shift_x = 0.172 * (ideal_beam_x - obs_beam_x)
        beam_shift_y = 0.172 * (ideal_beam_y - obs_beam_y)

        distance = (
            float(self._cif_header_dictionary["Detector_distance"].split()[0]) * 1000.0
        )

        wavelength = float(self._cif_header_dictionary["Wavelength"].split()[0])

        thickness = float(self._cif_header_dictionary["Silicon"].split()[2]) * 1000.0

        off_x = 184.9

        detector = Detector()
        root = detector.hierarchy()
        root.set_frame(
            x.elems,
            y.elems,
            (-distance * z + (beam_shift_x * x) + (beam_shift_y * y)).elems,
        )

        table = attenuation_coefficient.get_table("Si")
        mu = table.mu_at_angstrom(wavelength) / 10.0
        t0 = thickness
        px_mm = ParallaxCorrectedPxMmStrategy(mu, t0)

        self.coords = {}

        for j in range(24):
            shift_y = 195 + 17
            ymin, ymax = j * shift_y, j * shift_y + 195

            angle = math.pi * (-12.2 + 0.5 * 7.903 + j * (7.903 + 0.441)) / 180.0
            fast = matrix.col((1, 0, 0))
            slow = matrix.col((0, math.sin(angle), math.cos(angle)))
            normal = fast.cross(slow)

            row_origin = 250.0 * normal - off_x * fast - 16.8 * slow

            if not self._multi_panel:
                xmin, xmax = 0, 2463

                # OK two calls to add_panel here for detector like things => two
                # copies of the panel then? https://github.com/dials/dials/issues/189
                # ... this is also not the source of the leak

                # OBS! you need to set the panel to a root before set local frame...
                p = root.add_panel()
                p.set_type("SENSOR_PAD")
                p.set_name("row-%02d" % j)
                p.set_raw_image_offset((xmin, ymin))
                p.set_image_size((2463, 195))
                p.set_trusted_range((-1, 1000000))
                p.set_pixel_size((0.172, 0.172))
                p.set_local_frame(fast.elems, slow.elems, row_origin.elems)
                p.set_thickness(thickness)
                p.set_material("Si")
                p.set_mu(mu)
                p.set_px_mm_strategy(px_mm)
                p.set_raw_image_offset((xmin, ymin))
                self.coords[p.get_name()] = (xmin, ymin, xmax, ymax)

            else:
                shift_x = 487 + 7

                for i in range(5):
                    xmin, xmax = i * shift_x, i * shift_x + 487
                    origin = row_origin + i * (487 + 7) * 0.172 * fast

                    # OBS! you need to set the panel to a root before set local frame...
                    p = root.add_panel()
                    p.set_type("SENSOR_PAD")
                    p.set_name("row-%02d-col-%02d" % (j, i))
                    p.set_raw_image_offset((xmin, ymin))
                    p.set_image_size((487, 195))
                    p.set_trusted_range((-1, 1000000))
                    p.set_pixel_size((0.172, 0.172))
                    p.set_local_frame(fast.elems, slow.elems, origin.elems)
                    p.set_thickness(thickness)
                    p.set_material("Si")
                    p.set_mu(mu)
                    p.set_px_mm_strategy(px_mm)
                    p.set_raw_image_offset((xmin, ymin))
                    self.coords[p.get_name()] = (xmin, ymin, xmax, ymax)

        detector = self._mask_bad_modules(detector)
        return detector

    def _mask_bad_modules(self, detector):
        # Mask out known bad modules
        timestamp = get_pilatus_timestamp(self._cif_header_dictionary["timestamp"])
        # https://photon-science.desy.de/research/technical_groups/detectors/e190302/e199162/CharacterizationandCalibrationofPilatus.pdf
        nx = 487  # module pixels x
        ny = 195  # module pixels y
        cx = 60  # chip pixels x
        cy = 97  # chip pixels y
        dx = 7  # module gap size

        if timestamp > calendar.timegm((2020, 9, 8, 0, 0, 0)):
            # 2020 run 4
            # Detector serviced by Dectris, no bad modules
            pass
        elif timestamp > calendar.timegm((2020, 2, 21, 0, 0, 0)):
            # 2020 run 1
            # module @ row 1 column 4
            # blank chip @ row 15 column 4 (chip row 0, column 1)
            # module @ row 17 column 0
            # module @ row 17 column 4
            if self._multi_panel:
                detector[5 * 1 + 4].add_mask(0, 0, nx, ny)
                detector[5 * 17].add_mask(0, 0, nx, ny)
                detector[5 * 17 + 4].add_mask(0, 0, nx, ny)
                detector[5 * 15 + 4].add_mask(cx, 0, cx * 2, cy)
            else:
                detector[1].add_mask((nx + dx) * 4, 0, (nx + dx) * 4 + nx, ny)
                detector[17].add_mask(0, 0, nx, ny)
                detector[17].add_mask((nx + dx) * 4, 0, (nx + dx) * 4 + nx, ny)
                detector[15].add_mask((nx + dx) * 4 + cx, 0, (nx + dx) * 4 + cx * 2, cy)
        elif timestamp > calendar.timegm((2019, 11, 26, 0, 0, 0)):
            # 2019 run 5
            # module @ row 17 column 0
            # module @ row 17 column 4
            if self._multi_panel:
                detector[5 * 17].add_mask(0, 0, nx, ny)
                detector[5 * 17 + 4].add_mask(0, 0, nx, ny)
            else:
                detector[17].add_mask(0, 0, nx, ny)
                detector[17].add_mask((nx + dx) * 4, 0, (nx + dx) * 4 + nx, ny)
        elif timestamp > calendar.timegm((2019, 9, 3, 0, 0, 0)):
            # 2019 run 4
            # module @ row 15 column 2
            # module @ row 17 column 0
            if self._multi_panel:
                detector[5 * 15 + 2].add_mask(0, 0, nx, ny)
                detector[5 * 17].add_mask(0, 0, nx, ny)
            else:
                detector[15].add_mask((nx + dx) * 2, 0, (nx + dx) * 2 + nx, ny)
                detector[17].add_mask(0, 0, nx, ny)

        return detector

    def _read_cbf_image(self):
        start_tag = binascii.unhexlify("0c1a04d5")

        with self.open_file(self._image_file, "rb") as fh:
            data = fh.read()
        data_offset = data.find(start_tag) + 4
        cbf_header = self._parse_cbf_header(
            data[: data_offset - 4].decode("ascii", "ignore")
        )

        pixel_values = uncompress(
            packed=data[data_offset : data_offset + cbf_header["size"]],
            fast=cbf_header["fast"],
            slow=cbf_header["slow"],
        )

        return pixel_values

    def get_raw_data(self):
        if self._raw_data is None:
            raw_data = self._read_cbf_image()
            self._raw_data = []

            for panel in self.get_detector():
                xmin, ymin = panel.get_raw_image_offset()
                xmax = xmin + panel.get_image_size()[0]
                ymax = ymin + panel.get_image_size()[1]
                self._raw_data.append(raw_data[ymin:ymax, xmin:xmax])

        return tuple(self._raw_data)

    def get_goniometer_shadow_masker(self, goniometer=None):
        if goniometer is None:
            goniometer = self.get_goniometer()
        return GoniometerMaskerFactory.dls_i23_kappa(goniometer)

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header."""

        timestamp = get_pilatus_timestamp(self._cif_header_dictionary["timestamp"])
        # Goniometer changed from reverse phi to conventional rotation direction
        # on this date:
        # calendar.timegm(time.strptime('2016-04-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
        if timestamp < 1459468800:
            return self._goniometer_factory.single_axis_reverse()

        alpha = 50.0
        if "Phi" in self._cif_header_dictionary:
            phi_value = float(self._cif_header_dictionary["Phi"].split()[0])
        else:
            phi_value = 0.0

        if "Kappa" in self._cif_header_dictionary:
            kappa_value = float(self._cif_header_dictionary["Kappa"].split()[0])
        else:
            kappa_value = 0.0

        if "Omega" in self._cif_header_dictionary:
            omega_value = float(self._cif_header_dictionary["Omega"].split()[0])
        else:
            omega_value = 0.0

        return self._goniometer_factory.make_kappa_goniometer(
            alpha, omega_value, kappa_value, phi_value, "-y", "omega"
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDLS12M.understand(arg))
