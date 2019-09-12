"""
An implementation of the CBF image reader for Pilatus images, from the Pilatus
6M SN 100 currently on Diamond I04.
"""

from __future__ import absolute_import, division, print_function

import binascii
import sys

import libtbx
from cbflib_adaptbx import uncompress
from cctbx.eltbx import attenuation_coefficient
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask
from dxtbx.masking import GoniometerMaskerFactory
from dxtbx.model import ParallaxCorrectedPxMmStrategy
from dxtbx.model.detector import Detector

# Module positional offsets in x, y, in pixels - for the moment ignoring the
# rotational offsets as these are not well defined. To be honest these
# positional offsets are also not well defined as I do not know how they
# should be applied...

MODULE_OFFSETS_X = {
    (0, 0): -0.477546,
    (0, 1): 0.130578,
    (0, 2): 0.045041,
    (0, 3): -0.439872,
    (0, 4): -0.382077,
    (1, 0): 0.087405,
    (1, 1): 0.249597,
    (1, 2): 0.184265,
    (1, 3): 0.158342,
    (1, 4): 0.025225,
    (2, 0): -0.179892,
    (2, 1): -0.010974,
    (2, 2): -0.139207,
    (2, 3): 0.282851,
    (2, 4): -0.442219,
    (3, 0): -0.185027,
    (3, 1): 0.218601,
    (3, 2): 0.092585,
    (3, 3): 0.35862,
    (3, 4): -0.29161,
    (4, 0): 0.145368,
    (4, 1): 0.609289,
    (4, 2): 0.396265,
    (4, 3): 0.41625,
    (4, 4): 0.07152,
    (5, 0): 0.247142,
    (5, 1): 0.046563,
    (5, 2): 0.248714,
    (5, 3): -0.044628,
    (5, 4): -0.391509,
    (6, 0): 0.516643,
    (6, 1): 0.358453,
    (6, 2): 0.069219,
    (6, 3): 0.095861,
    (6, 4): -0.167403,
    (7, 0): -0.381352,
    (7, 1): -0.35338,
    (7, 2): 0.348656,
    (7, 3): 0.024543,
    (7, 4): 0.328706,
    (8, 0): 0.150886,
    (8, 1): 0.244987,
    (8, 2): -0.102911,
    (8, 3): 0.16633,
    (8, 4): 0.386622,
    (9, 0): 0.037924,
    (9, 1): 0.314392,
    (9, 2): 0.238818,
    (9, 3): 0.815028,
    (9, 4): -0.048818,
    (10, 0): -0.670524,
    (10, 1): -0.304119,
    (10, 2): 0.252284,
    (10, 3): -0.05485,
    (10, 4): -0.355264,
    (11, 0): -0.404947,
    (11, 1): -0.020622,
    (11, 2): 0.648473,
    (11, 3): -0.277175,
    (11, 4): -0.711951,
}

MODULE_OFFSETS_Y = {
    (0, 0): -0.494797,
    (0, 1): -0.212976,
    (0, 2): 0.085351,
    (0, 3): 0.35494,
    (0, 4): 0.571189,
    (1, 0): -0.421708,
    (1, 1): 0.061914,
    (1, 2): 0.238996,
    (1, 3): 0.146692,
    (1, 4): 0.407145,
    (2, 0): -0.313212,
    (2, 1): -0.225025,
    (2, 2): 0.031613,
    (2, 3): -0.047839,
    (2, 4): 0.42716,
    (3, 0): -0.361193,
    (3, 1): 0.057663,
    (3, 2): 0.022357,
    (3, 3): 0.062717,
    (3, 4): 0.150611,
    (4, 0): 0.035511,
    (4, 1): -0.271567,
    (4, 2): 0.007761,
    (4, 3): -0.124021,
    (4, 4): 0.093017,
    (5, 0): -0.238897,
    (5, 1): -0.179724,
    (5, 2): -0.113608,
    (5, 3): 0.017841,
    (5, 4): -0.012933,
    (6, 0): -0.166337,
    (6, 1): -0.272922,
    (6, 2): -0.194665,
    (6, 3): -0.058535,
    (6, 4): -0.405404,
    (7, 0): -0.318824,
    (7, 1): -0.311276,
    (7, 2): -0.205223,
    (7, 3): -0.292664,
    (7, 4): -0.474762,
    (8, 0): -0.039504,
    (8, 1): -0.239887,
    (8, 2): -0.343485,
    (8, 3): -0.459429,
    (8, 4): -0.426901,
    (9, 0): -0.187805,
    (9, 1): 0.282727,
    (9, 2): -0.601164,
    (9, 3): -0.467605,
    (9, 4): -0.589271,
    (10, 0): 0.028311,
    (10, 1): -0.391571,
    (10, 2): -0.463112,
    (10, 3): -0.358092,
    (10, 4): -0.285396,
    (11, 0): 0.01863,
    (11, 1): -0.380099,
    (11, 2): -0.234953,
    (11, 3): -0.593992,
    (11, 4): -0.801247,
}


class FormatCBFMiniPilatusDLS6MSN100(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 100 @ DLS."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 60-0100 Diamond" in header
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

        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)
        self._multi_panel = kwargs.get("multi_panel", False)
        super(FormatCBFMiniPilatusDLS6MSN100, self).__init__(image_file, **kwargs)

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header, though for miniCBF
        there are limited options for this."""

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

        phi = (1.0, 0.0, 0.0)
        kappa = (0.914, 0.279, -0.297)
        omega = (1.0, 0.0, 0.0)
        axes = flex.vec3_double((phi, kappa, omega))
        angles = flex.double((phi_value, kappa_value, omega_value))
        names = flex.std_string(("GON_PHI", "GON_KAPPA", "GON_OMEGA"))
        return self._goniometer_factory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    def _detector(self):
        """Detector model, allowing for small offsets in the positions of 60
        detector modules"""

        distance = float(self._cif_header_dictionary["Detector_distance"].split()[0])

        beam_xy = (
            self._cif_header_dictionary["Beam_xy"]
            .replace("(", "")
            .replace(")", "")
            .replace(",", "")
            .split()[:2]
        )

        beam_x, beam_y = map(float, beam_xy)

        wavelength = float(self._cif_header_dictionary["Wavelength"].split()[0])

        pixel_xy = (
            self._cif_header_dictionary["Pixel_size"]
            .replace("m", "")
            .replace("x", "")
            .split()
        )

        pixel_x, pixel_y = map(float, pixel_xy)

        thickness = float(self._cif_header_dictionary["Silicon"].split()[2]) * 1000.0

        nx = int(self._cif_header_dictionary["X-Binary-Size-Fastest-Dimension"])
        ny = int(self._cif_header_dictionary["X-Binary-Size-Second-Dimension"])

        overload = int(self._cif_header_dictionary["Count_cutoff"].split()[0])
        underload = -1

        # take into consideration here the thickness of the sensor also the
        # wavelength of the radiation (which we have in the same file...)
        table = attenuation_coefficient.get_table("Si")
        mu = table.mu_at_angstrom(wavelength) / 10.0
        t0 = thickness

        # FIXME would also be very nice to be able to take into account the
        # misalignment of the individual modules given the calibration...

        # single detector or multi-module detector

        pixel_x *= 1000.0
        pixel_y *= 1000.0
        distance *= 1000.0

        if not self._multi_panel:
            detector = self._detector_factory.simple(
                "PAD",
                distance,
                (beam_x * pixel_x, beam_y * pixel_y),
                "+x",
                "-y",
                (pixel_x, pixel_y),
                (nx, ny),
                (underload, overload),
                [],
                ParallaxCorrectedPxMmStrategy(mu, t0),
            )

            for f0, f1, s0, s1 in determine_pilatus_mask(detector):
                detector[0].add_mask(f0 - 1, s0 - 1, f1, s1)

            detector[0].set_thickness(thickness)
            detector[0].set_material("Si")
            detector[0].set_mu(mu)

            return detector

        # got to here means 60-panel version
        d = Detector()

        beam_centre = matrix.col((beam_x * pixel_x, beam_y * pixel_y, 0))

        fast = matrix.col((1.0, 0.0, 0.0))
        slow = matrix.col((0.0, -1.0, 0.0))
        s0 = matrix.col((0, 0, -1))
        origin = (distance * s0) - (fast * beam_centre[0]) - (slow * beam_centre[1])

        root = d.hierarchy()
        root.set_local_frame(fast.elems, slow.elems, origin.elems)

        xmins = [0, 494, 988, 1482, 1976]
        xmaxes = [487, 981, 1475, 1969, 2463]
        ymins = [0, 212, 424, 636, 848, 1060, 1272, 1484, 1696, 1908, 2120, 2332]
        ymaxes = [195, 407, 619, 831, 1043, 1255, 1467, 1679, 1891, 2103, 2315, 2527]

        self.coords = {}

        fast = matrix.col((1.0, 0.0, 0.0))
        slow = matrix.col((0.0, 1.0, 0.0))
        panel_idx = 0
        for ymin, ymax in zip(ymins, ymaxes):
            for xmin, xmax in zip(xmins, xmaxes):
                xmin_mm = xmin * pixel_x
                ymin_mm = ymin * pixel_y

                origin_panel = fast * xmin_mm + slow * ymin_mm

                panel_name = "Panel%d" % panel_idx
                panel_idx += 1

                p = d.add_panel()
                p.set_type("SENSOR_PAD")
                p.set_name(panel_name)
                p.set_raw_image_offset((xmin, ymin))
                p.set_image_size((xmax - xmin, ymax - ymin))
                p.set_trusted_range((underload, overload))
                p.set_pixel_size((pixel_x, pixel_y))
                p.set_thickness(thickness)
                p.set_material("Si")
                p.set_mu(mu)
                p.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
                p.set_local_frame(fast.elems, slow.elems, origin_panel.elems)
                p.set_raw_image_offset((xmin, ymin))
                self.coords[panel_name] = (xmin, ymin, xmax, ymax)

        return d

    def get_goniometer_shadow_masker(self, goniometer=None):
        return GoniometerMaskerFactory.mini_kappa(goniometer)

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
        if not self._multi_panel:
            return super(FormatCBFMiniPilatusDLS6MSN100, self).get_raw_data()

        if self._raw_data is None:
            raw_data = self._read_cbf_image()
            self._raw_data = []
            d = self.get_detector()
            for panel in d:
                xmin, ymin, xmax, ymax = self.coords[panel.get_name()]
                self._raw_data.append(raw_data[ymin:ymax, xmin:xmax])
            self._raw_data = tuple(self._raw_data)
        return self._raw_data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDLS6MSN100.understand(arg))
