"""
An implementation of the CBF image reader for ADSC images, from the ADSC
HF-4M SN H401 currently on APS sector 24 (NE-CAT).
Located in dxtbx/format
"""


import binascii
import calendar
import re
import sys
import time

from cctbx.eltbx import attenuation_coefficient
from iotbx.detectors.adsc_minicbf import ADSCHF4MImage

from dxtbx.ext import uncompress
from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.model import ParallaxCorrectedPxMmStrategy


def get_adsc_timestamp(timestamp):
    # reduce duplicate to single underscore
    timestamp = re.sub("_+", "_", timestamp)

    for format in ["%a_%b_%d_%H:%M:%S_%Y"]:

        try:
            struct_time = time.strptime(timestamp, format)
            return calendar.timegm(struct_time)

        except Exception:
            pass

    raise RuntimeError("timestamp %s not recognised" % timestamp)


class FormatCBFMiniADSCHF4M(FormatCBFMini):
    """A class for reading mini CBF format ADSC images for HF-4M @ NE-CAT."""

    @staticmethod
    def understand(image_file):

        header = FormatCBFMini.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "# Detector" in record and "ADSC" in record and "HF-4M" in header:
                return True
            if "_array_data.header_convention" in record and "ADSC" in record:
                return True

        return False

    def _detector(self):
        """Return a model for a simple detector, presuming no one has
        one of these on a two-theta stage. Assert that the beam centre is
        provided in the Mosflm coordinate frame."""

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
        # misalignment of the individual modules given the calibration information...

        detector = self._detector_factory.simple(
            "PAD",
            distance * 1000.0,
            (beam_x * pixel_x * 1000.0, beam_y * pixel_y * 1000.0),
            "+x",
            "-y",
            (1000 * pixel_x, 1000 * pixel_y),
            (nx, ny),
            (underload, overload),
            [],
            ParallaxCorrectedPxMmStrategy(mu, t0),
        )
        """
    for f0, f1, s0, s1 in determine_pilatus_mask(detector):
      detector[0].add_mask(f0-1, s0-1, f1, s1)
    """
        return detector

    def _beam(self):
        """Return a simple model for the beam."""

        wavelength = float(self._cif_header_dictionary["Wavelength"].split()[0])

        beam = self._beam_factory.simple(wavelength)

        try:
            flux = float(self._cif_header_dictionary["Flux"].split()[0])
            beam.set_flux(flux)
        except KeyError:
            pass

        try:
            transmission = float(self._cif_header_dictionary["Transmission"].split()[0])
            beam.set_transmission(transmission)
        except KeyError:
            pass

        return beam

    def _scan(self):
        """Return the scan information for this image."""

        exposure_time = float(self._cif_header_dictionary["Exposure_period"].split()[0])

        osc_start = float(self._cif_header_dictionary["Start_angle"].split()[0])
        osc_range = float(self._cif_header_dictionary["Angle_increment"].split()[0])

        timestamp = get_adsc_timestamp(self._cif_header_dictionary["timestamp"])

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, timestamp
        )

    def detectorbase_start(self):

        self.detectorbase = ADSCHF4MImage(self._image_file)
        self.detectorbase.readHeader()

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
            data = self._read_cbf_image()
            self._raw_data = data

        return self._raw_data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniADSCHF4M.understand(arg))
