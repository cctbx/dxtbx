"""An implementation of the CBF image reader for Rigaku Hypix 6000 images"""


import time

from cctbx.eltbx import attenuation_coefficient
from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
from scitbx import matrix

from dials.array_family import flex

from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.model import ParallaxCorrectedPxMmStrategy


class FormatCBFMiniRigakuHypix6000(FormatCBFMini):
    """Support for mini CBF format RigakuHypix6000 images."""

    def __init__(self, image_file, **kwargs):
        super().__init__(image_file, **kwargs)

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMini.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "# Detector" in record and "HyPix-6000" in record:
                return True

        return False

    def _start(self):
        """Open the image file, read the image header, copy it into a
        dictionary for future reference."""

        super()._start()
        cif_header = FormatCBFMini.get_cbf_header(self._image_file)

        self._cif_header_dictionary = {}

        for record in cif_header.split("\n"):
            if record[:1] != "#":
                continue

            if len(record[1:].split()) <= 2 and record.count(":") == 2:
                self._cif_header_dictionary["timestamp"] = record[1:].strip()
                continue

            tokens = record.replace('"', "").split()
            self._cif_header_dictionary[tokens[1]] = " ".join(tokens[2:])

        for record in self._mime_header.split("\n"):
            if not record.strip():
                continue
            token, value = record.split(":")
            self._cif_header_dictionary[token.strip()] = value.strip()

        _X = tuple(
            map(float, self._cif_header_dictionary["Rotation_axis_vector"].split())
        )
        _Z = tuple(
            map(float, self._cif_header_dictionary["Incident_beam_vector"].split())
        )

        self._R = align_reference_frame(_X, (1, 0, 0), _Z, (0, 0, 1))
        self._X = self._R * _X
        self._Z = self._R * _Z

        self._multi_panel = False

    def _detector(self):
        distance = float(self._cif_header_dictionary["Detector_distance"].split()[0])

        beam_xy = (
            self._cif_header_dictionary["Beam_xy"]
            .replace("(", "")
            .replace(")", "")
            .replace(",", "")
            .split()[:2]
        )

        wavelength = float(self._cif_header_dictionary["Wavelength"].split()[0])

        beam_x, beam_y = map(float, beam_xy)

        pixel_xy = (
            self._cif_header_dictionary["Pixel_size"]
            .replace("m", "")
            .replace("x", "")
            .split()
        )

        pixel_x, pixel_y = map(float, pixel_xy)

        if "Silicon" in self._cif_header_dictionary:
            thickness = (
                float(self._cif_header_dictionary["Silicon"].split()[2]) * 1000.0
            )
            material = "Si"
            sensor = "PAD"
        else:
            raise RuntimeError("Unsupported detector material: please contact authors")

        nx = int(self._cif_header_dictionary["X-Binary-Size-Fastest-Dimension"])
        ny = int(self._cif_header_dictionary["X-Binary-Size-Second-Dimension"])

        # no idea what the proper limit is here
        overload = 0xFFFFFF
        underload = -1

        table = attenuation_coefficient.get_table(material)
        mu = table.mu_at_angstrom(wavelength) / 10.0
        t0 = thickness
        px_mm = ParallaxCorrectedPxMmStrategy(mu, t0)

        # yes, fast and slow are inverted
        fast = self._R * tuple(
            map(float, self._cif_header_dictionary["Detector_slow_axis_vector"].split())
        )
        slow = self._R * tuple(
            map(float, self._cif_header_dictionary["Detector_fast_axis_vector"].split())
        )

        # fast and slow have already been rotated but origin needs rotating
        origin = (
            -(beam_x * pixel_x * 1000.0 * fast)
            - (beam_y * pixel_y * 1000.0 * slow)
            - (distance * 1000 * matrix.col((0, 0, 1)))
        ).rotate_around_origin(
            self._X,
            float(self._cif_header_dictionary["Detector_2theta"].split()[0]),
            deg=True,
        )

        detector = self._detector_factory.make_detector(
            sensor,
            fast,
            slow,
            origin,
            (1000 * pixel_x, 1000 * pixel_y),
            (nx, ny),
            (underload, overload),
            px_mm=px_mm,
            mu=mu,
        )

        detector[0].set_thickness(thickness)
        detector[0].set_material(material)
        detector[0].set_mu(mu)

        return detector

    def _scan(self):
        exposure_time = float(self._cif_header_dictionary["TIME"].split()[0])

        osc_start = float(self._cif_header_dictionary["Start_angle"].split()[0])
        osc_range = float(self._cif_header_dictionary["Angle_increment"].split()[0])

        timestamp = time.mktime(time.strptime(self._cif_header_dictionary["TIMESTAMP"]))

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, timestamp
        )

    def _goniometer(self):
        phi = self._R * tuple(
            map(float, self._cif_header_dictionary["Phi_axis_vector"].split())
        )
        kappa = self._R * tuple(
            map(float, self._cif_header_dictionary["Kappa_axis_vector"].split())
        )
        omega = self._R * tuple(
            map(float, self._cif_header_dictionary["Omega_axis_vector"].split())
        )

        axes = flex.vec3_double((phi, kappa, omega))

        _phi = float(self._cif_header_dictionary["Phi"].split()[0])
        _kappa = float(self._cif_header_dictionary["Kappa"].split()[0])
        _omega = float(self._cif_header_dictionary["Omega"].split()[0])

        angles = flex.double((_phi, _kappa, _omega))

        names = flex.std_string(("PHI", "KAPPA", "OMEGA"))

        return self._goniometer_factory.multi_axis(axes, angles, names, 2)

    def _beam(self):
        wavelength = float(self._cif_header_dictionary["Wavelength"].split()[0])

        beam = self._beam_factory.make_beam(
            sample_to_source=self._Z, wavelength=wavelength
        )

        return beam

    def get_vendortype(self):
        return "Rigaku"
