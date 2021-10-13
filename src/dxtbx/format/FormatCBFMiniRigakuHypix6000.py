"""An implementation of the CBF image reader for Rigaku Hypix 6000 images"""


import time

from cctbx.eltbx import attenuation_coefficient

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

        fast = 1, 0, 0
        slow = 0, -1, 0
        origin = -beam_x * pixel_x * 1000.0, beam_y * pixel_y * 1000.0, -distance * 1000

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

    def get_vendortype(self):
        return "Rigaku"

    @staticmethod
    def as_file(
        detector,
        beam,
        gonio,
        scan,
        data,
        path,
        header_convention="PILATUS_1.2",
        det_type="PILATUS3 6M",
    ):
        FormatCBFMini.as_file(
            detector, beam, gonio, scan, data, path, header_convention, det_type
        )

    def get_raw_data(self):
        if not self._multi_panel:
            return super().get_raw_data()

        if self._raw_data is None:
            raw_data = self._read_cbf_image()
            self._raw_data = []
            d = self.get_detector()
            for panel in d:
                xmin, ymin, xmax, ymax = self.coords[panel.get_name()]
                self._raw_data.append(raw_data[ymin:ymax, xmin:xmax])
            self._raw_data = tuple(self._raw_data)
        return self._raw_data
