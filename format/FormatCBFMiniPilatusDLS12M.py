#!/usr/bin/env python
# FormatCBFMiniPilatusDLS12M.py
#
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images, for the P12M-DLS

from __future__ import division
from __future__ import print_function

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.model import ParallaxCorrectedPxMmStrategy

__mask = None

# if group_rows == True, then interpret data as 24 panels, where each row
# of 5 panels is grouped as one "panel"
# elif group_rows == False, then interpret data as 120 panels,
# 24 rows * 5 columns - bodge through ENV variable in 1st cut...

import os

if "P12M_120_PANEL" in os.environ:
    group_rows = False
else:
    group_rows = True

if "DXTBX_ENABLE_SHADOWING" in os.environ:
    enable_shadowing = True
else:
    enable_shadowing = False


def read_mask():
    global __mask
    if not __mask:
        import os
        import cPickle as pickle

        source_dir = os.path.split(__file__)[0]
        mask_file = os.path.join(source_dir, "FormatCBFMiniPilatusDLS12M.pickle")
        __mask = pickle.load(open(mask_file, "rb"))
    return __mask


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

    def __init__(self, image_file):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        assert self.understand(image_file)

        FormatCBFMiniPilatus.__init__(self, image_file)

        self._raw_data = None

        return

    def _detector(self):

        # module positions from detector blueprints - modelling at the moment as
        # 24 modules, each consisting of 5 sensors (the latter is ignored)

        from dxtbx.model import Detector
        from scitbx import matrix
        import math

        x = matrix.col((-1, 0, 0))
        y = matrix.col((0, 1, 0))
        z = matrix.col((0, 0, 1))

        obs_beam_y = 2587
        ideal_beam_y = 2594
        beam_shift_y = 0.172 * (2594 - 2587)

        distance = (
            float(self._cif_header_dictionary["Detector_distance"].split()[0]) * 1000.0
        )

        wavelength = float(self._cif_header_dictionary["Wavelength"].split()[0])

        thickness = float(self._cif_header_dictionary["Silicon"].split()[2]) * 1000.0

        # for longer wavelength data sets move 192.3 below to 184.9
        if wavelength < 1.128:
            off_x = 191.9
        else:
            off_x = 184.9

        z += beam_shift_y * y

        detector = Detector()
        root = detector.hierarchy()
        root.set_frame(x.elems, y.elems, (-distance * z).elems)

        from cctbx.eltbx import attenuation_coefficient

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

            if group_rows:
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

        return detector

    def read_cbf_image(self, cbf_image):
        from cbflib_adaptbx import uncompress
        import binascii

        start_tag = binascii.unhexlify("0c1a04d5")

        data = self.open_file(cbf_image, "rb").read()
        data_offset = data.find(start_tag) + 4
        cbf_header = data[: data_offset - 4]

        fast = 0
        slow = 0
        length = 0

        for record in cbf_header.split("\n"):
            if "X-Binary-Size-Fastest-Dimension" in record:
                fast = int(record.split()[-1])
            elif "X-Binary-Size-Second-Dimension" in record:
                slow = int(record.split()[-1])
            elif "X-Binary-Number-of-Elements" in record:
                length = int(record.split()[-1])
            elif "X-Binary-Size:" in record:
                size = int(record.split()[-1])

        assert length == fast * slow

        pixel_values = uncompress(
            packed=data[data_offset : data_offset + size], fast=fast, slow=slow
        )

        isel = read_mask()
        pixel_values.as_1d().set_selected(isel, -2)

        return pixel_values

    def get_raw_data(self):
        if self._raw_data is None:
            raw_data = self.read_cbf_image(self._image_file)
            self._raw_data = []

            for panel in self.get_detector():
                xmin, ymin = panel.get_raw_image_offset()
                xmax = xmin + panel.get_image_size()[0]
                ymax = ymin + panel.get_image_size()[1]
                self._raw_data.append(raw_data[ymin:ymax, xmin:xmax])

        return tuple(self._raw_data)

    if enable_shadowing:

        def get_goniometer_shadow_masker(self):
            from dials.util.masking import GoniometerShadowMaskGenerator
            from scitbx.array_family import flex
            import math

            coords = flex.vec3_double(((0, 0, 0),))

            r = flex.double(19, 40)
            alpha = flex.double_range(0, 190, step=10) * math.pi / 180
            x = flex.double(r.size(), 107.61)
            y = -r * flex.sin(alpha)
            z = -r * flex.cos(alpha)
            coords.extend(flex.vec3_double(x, y, z))

            coords.extend(
                flex.vec3_double(
                    (
                        # fixed
                        (107.49, 7.84, 39.49),
                        (107.39, 15.69, 38.97),
                        (107.27, 23.53, 38.46),
                        (107.16, 31.37, 37.94),
                        (101.76, 33.99, 36.25),
                        (96.37, 36.63, 34.56),
                        (90.98, 39.25, 33.00),
                        (85.58, 41.88, 31.18),
                        (80.89, 47.06, 31.00),
                        (76.55, 51.51, 31.03),
                        (72.90, 55.04, 31.18),
                        (66.86, 60.46, 31.67),
                        (62.10, 64.41, 32.25),
                    )
                )
            )

            r = flex.double(19, 33)
            alpha = flex.double_range(180, 370, step=10) * math.pi / 180
            x = flex.sqrt(flex.pow2(r * flex.sin(alpha)) + 89.02 ** 2) * flex.cos(
                (50 * math.pi / 180) - flex.atan(r / 89.02 * flex.sin(alpha))
            )
            y = flex.sqrt(flex.pow2(r * flex.sin(alpha)) + 89.02 ** 2) * flex.sin(
                (50 * math.pi / 180) - flex.atan(r / 89.02 * flex.sin(alpha))
            )
            z = -r * flex.cos(alpha)
            coords.extend(flex.vec3_double(x, y, z))

            coords.extend(
                flex.vec3_double(
                    (
                        # fixed
                        (62.10, 64.41, -32.25),
                        (66.86, 60.46, -31.67),
                        (72.90, 55.04, -31.18),
                        (76.55, 51.51, -31.03),
                        (80.89, 47.06, -31.00),
                        (85.58, 41.88, -31.18),
                        (90.98, 39.25, -33.00),
                        (96.37, 36.63, -34.56),
                        (101.76, 33.99, -36.25),
                        (107.16, 31.37, -37.94),
                        (107.27, 23.53, -38.46),
                        (107.39, 15.69, -38.97),
                        (107.49, 7.84, -39.49),
                        (107.61, 0.00, -40.00),
                    )
                )
            )

            gonio = self.get_goniometer()
            return GoniometerShadowMaskGenerator(
                gonio, coords, flex.size_t(len(coords), 0)
            )

        def get_mask(self):
            gonio_masker = self.get_goniometer_shadow_masker()
            scan = self.get_scan()
            detector = self.get_detector()
            mask = super(FormatCBFMiniPilatusDLS6MSN100, self).get_mask()
            shadow_mask = gonio_masker.get_mask(detector, scan.get_oscillation()[0])
            assert len(mask) == len(shadow_mask)
            for m, sm in zip(mask, shadow_mask):
                m &= sm

            return mask

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header."""

        from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp

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
            alpha, omega_value, kappa_value, phi_value, "+y", "omega"
        )


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDLS12M.understand(arg))
