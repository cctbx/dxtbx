"""An implementation of the CBF image reader for Eiger images"""


import binascii
import os
import sys

from cctbx.eltbx import attenuation_coefficient
from iotbx.detectors.eiger_minicbf import EigerCBFImage

from dxtbx.ext import uncompress
from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp
from dxtbx.format.FormatPilatusHelpers import determine_eiger_mask
from dxtbx.format.FormatPilatusHelpers import get_vendortype_eiger as gv
from dxtbx.model import ParallaxCorrectedPxMmStrategy

if "DXTBX_OVERLOAD_SCALE" in os.environ:
    dxtbx_overload_scale = float(os.environ["DXTBX_OVERLOAD_SCALE"])
else:
    dxtbx_overload_scale = 1


class FormatCBFMiniEiger(FormatCBFMini):
    """A class for reading mini CBF format Eiger images, and correctly
    constructing a model for the experiment from this."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Eiger mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMini.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "# detector" in record.lower() and "eiger" in record.lower():
                return True

        return False

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
        elif "CdTe" in self._cif_header_dictionary:
            thickness = float(self._cif_header_dictionary["CdTe"].split()[2]) * 1000.0
            material = "CdTe"
        else:
            thickness = 0.450
            material = "Si"

        nx = int(self._cif_header_dictionary["X-Binary-Size-Fastest-Dimension"])
        ny = int(self._cif_header_dictionary["X-Binary-Size-Second-Dimension"])

        if "Count_cutoff" in self._cif_header_dictionary:
            overload = int(self._cif_header_dictionary["Count_cutoff"].split()[0])
        else:
            # missing from data transformed with GPhL converter - dials#376
            overload = 100000000
        underload = -1

        try:
            identifier = self._cif_header_dictionary["Detector"].encode()
        except KeyError:
            identifier = "Unknown Eiger"

        table = attenuation_coefficient.get_table(material)
        mu = table.mu_at_angstrom(wavelength) / 10.0
        t0 = thickness

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
            px_mm=ParallaxCorrectedPxMmStrategy(mu, t0),
            mu=mu,
        )

        for f0, f1, s0, s1 in determine_eiger_mask(detector):
            detector[0].add_mask(f0 - 1, s0 - 1, f1, s1)

        for panel in detector:
            panel.set_thickness(thickness)
            panel.set_material(material)
            panel.set_identifier(identifier)
            panel.set_mu(mu)

        return detector

    def _beam(self):
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
        exposure_time = float(self._cif_header_dictionary["Exposure_period"].split()[0])
        osc_start = float(self._cif_header_dictionary["Start_angle"].split()[0])
        osc_range = float(self._cif_header_dictionary["Angle_increment"].split()[0])

        if "timestamp" in self._cif_header_dictionary:
            timestamp = get_pilatus_timestamp(self._cif_header_dictionary["timestamp"])
        else:
            timestamp = 0.0

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, timestamp
        )

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

    def detectorbase_start(self):
        self.detectorbase = EigerCBFImage(self._image_file)
        self.detectorbase.readHeader()

    def get_vendortype(self):
        return gv(self.get_detector())


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniEiger.understand(arg))
