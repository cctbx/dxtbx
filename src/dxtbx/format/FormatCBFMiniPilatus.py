"""An implementation of the CBF image reader for Pilatus images"""


import os
import sys

from cctbx.eltbx import attenuation_coefficient
from scitbx import matrix

from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp
from dxtbx.format.FormatPilatusHelpers import _DetectorDatabase, determine_pilatus_mask
from dxtbx.format.FormatPilatusHelpers import get_vendortype as gv
from dxtbx.model import Detector, ParallaxCorrectedPxMmStrategy


class FormatCBFMiniPilatus(FormatCBFMini):
    """A class for reading mini CBF format Pilatus images, and correctly
    constructing a model for the experiment from this."""

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        self._multi_panel = kwargs.get("multi_panel", False)
        super().__init__(image_file, **kwargs)

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        if "ENABLE_PHOTON_FACTORY_TWO_EIGER" in os.environ:
            return False

        header = FormatCBFMini.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "# Detector" in record:
                if "EIGER" in record.upper():
                    return False
                if "TIMEPIX" in record.upper():
                    return False

        for record in header.split("\n"):
            if "_array_data.header_convention" in record and "PILATUS" in record:
                return True
            if "_array_data.header_convention" in record and "SLS" in record:
                return True
            if "_array_data.header_convention" in record and "?" in record:
                return True
            if "# Detector" in record and "PILATUS" in record:  # CBFlib v0.8.0 allowed
                return True

        return False

    def _detector(self):
        """Return a model for a simple detector, presuming no one has
        one of these on a two-theta stage. Assert that the beam centre is
        provided in the Mosflm coordinate frame."""

        if not self._multi_panel:
            detector = FormatCBFMini._detector(self)
            for f0, f1, s0, s1 in determine_pilatus_mask(detector):
                detector[0].add_mask(f0 - 1, s0 - 1, f1, s1)
            return detector

        # got to here means 60-panel version
        d = Detector()

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

        beam_centre = matrix.col((beam_x * pixel_x, beam_y * pixel_y, 0))

        fast = matrix.col((1.0, 0.0, 0.0))
        slow = matrix.col((0.0, -1.0, 0.0))
        s0 = matrix.col((0, 0, -1))
        origin = (distance * s0) - (fast * beam_centre[0]) - (slow * beam_centre[1])

        root = d.hierarchy()
        root.set_local_frame(fast.elems, slow.elems, origin.elems)

        det = _DetectorDatabase["Pilatus"]

        # Edge dead areas not included, only gaps between modules matter
        n_fast, remainder = divmod(nx, det.module_size_fast)
        assert (n_fast - 1) * det.gap_fast == remainder

        n_slow, remainder = divmod(ny, det.module_size_slow)
        assert (n_slow - 1) * det.gap_slow == remainder

        mx = det.module_size_fast
        my = det.module_size_slow
        dx = det.gap_fast
        dy = det.gap_slow

        xmins = [(mx + dx) * i for i in range(n_fast)]
        xmaxes = [mx + (mx + dx) * i for i in range(n_fast)]
        ymins = [(my + dy) * i for i in range(n_slow)]
        ymaxes = [my + (my + dy) * i for i in range(n_slow)]

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

    def _scan(self):
        """Return the scan information for this image."""

        exposure_time = float(self._cif_header_dictionary["Exposure_period"].split()[0])

        osc_start = float(self._cif_header_dictionary["Start_angle"].split()[0])
        osc_range = float(self._cif_header_dictionary["Angle_increment"].split()[0])

        timestamp = get_pilatus_timestamp(self._cif_header_dictionary["timestamp"])

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, timestamp
        )

    def get_vendortype(self):
        return gv(self.get_detector())

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


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatus.understand(arg))
