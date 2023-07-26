"""
Base implementation of miniCBF format - as used with Dectris detectors -
this will read the header and populate a dictionary of the keyword / value
pairs.
"""


from __future__ import annotations

import binascii
import datetime
import os
import pathlib
import sys

import pycbf

from boost_adaptbx.boost.python import streambuf
from cctbx import factor_ev_angstrom
from cctbx.eltbx import attenuation_coefficient
from iotbx.detectors.pilatus_minicbf import PilatusImage
from scitbx.array_family import flex

from dxtbx.ext import read_int32, uncompress
from dxtbx.format.FormatCBF import FormatCBF
from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp
from dxtbx.format.FormatCBFMultiTile import cbf_wrapper
from dxtbx.model import ParallaxCorrectedPxMmStrategy, SimplePxMmStrategy

dxtbx_overload_scale = float(os.getenv("DXTBX_OVERLOAD_SCALE", "1"))


class FormatCBFMini(FormatCBF):
    """An image reading class for mini CBF format images i.e. those from
    Dectris, which will read the header into a dictionary."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        header = FormatCBF.get_cbf_header(image_file)

        if "_diffrn.id" in header and "_diffrn_source" in header:
            return False

        def one_of_these_in(record):
            these = [
                "PILATUS",
                "SLS",
                "SSRL",
                "?",
                "XDS special",
                "GENERIC_MINI",  # intended for simulated PAD data, non-Pilatus array size
            ]
            for convention in these:
                if convention in record:
                    return True
            return False

        for record in header.split("\n"):
            if "_array_data.header_convention" in record and one_of_these_in(record):
                return True
            if "# Detector" in record and "PILATUS" in record:  # CBFlib v0.8.0 allowed
                return True
            if "# Detector" in record and "ADSC" in record and "HF-4M" in header:
                return True

        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        self._raw_data = None
        super().__init__(image_file, **kwargs)

    @staticmethod
    def _get_timestamp_from_raw_header(
        header: str | list[str],
    ) -> datetime.datetime | None:
        """Given a raw header, or lines from, attempt to extract the timestamp field"""
        if isinstance(header, str):
            header = header.splitlines()
        timestamp = None
        for record in header:
            if len(record[1:].split()) <= 2 and record.count(":") == 2:
                timestamp = datetime.datetime.fromisoformat(record[1:].strip())
                break
        return timestamp

    def _start(self):
        """Open the image file, read the image header, copy it into a
        dictionary for future reference."""

        super()._start()
        cif_header = FormatCBF.get_cbf_header(self._image_file)

        self._cif_header_dictionary = {}

        for record in cif_header.split("\n"):
            if record[:1] != "#":
                continue

            if len(record[1:].split()) <= 2 and record.count(":") == 2:
                self._cif_header_dictionary["timestamp"] = record[1:].strip()
                continue

            tokens = record.replace("=", "").replace(":", "").split()
            self._cif_header_dictionary[tokens[1]] = " ".join(tokens[2:])

        for record in self._mime_header.split("\n"):
            if not record.strip():
                continue
            token, value = record.split(":")
            self._cif_header_dictionary[token.strip()] = value.strip()

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

        if "Silicon" in self._cif_header_dictionary:
            thickness = (
                float(self._cif_header_dictionary["Silicon"].split()[2]) * 1000.0
            )
            material = "Si"
            sensor = "PAD"
        elif "CdTe" in self._cif_header_dictionary:
            thickness = float(self._cif_header_dictionary["CdTe"].split()[2]) * 1000.0
            material = "CdTe"
            sensor = "PAD"
        elif "CCD" in self._cif_header_dictionary:
            thickness = 0
            material = None
            sensor = "CCD"
        else:
            thickness = 0
            material = None
            sensor = None

        nx = int(self._cif_header_dictionary["X-Binary-Size-Fastest-Dimension"])
        ny = int(self._cif_header_dictionary["X-Binary-Size-Second-Dimension"])

        overload = dxtbx_overload_scale * int(
            self._cif_header_dictionary["Count_cutoff"].split()[0]
        )
        minimum_trusted_value = 0

        if material is not None:
            # take into consideration here the thickness of the sensor also the
            # wavelength of the radiation (which we have in the same file...)
            table = attenuation_coefficient.get_table(material)
            mu = table.mu_at_angstrom(wavelength) / 10.0
            t0 = thickness
            px_mm = ParallaxCorrectedPxMmStrategy(mu, t0)
        else:
            px_mm = SimplePxMmStrategy()

        detector = self._detector_factory.simple(
            sensor,
            distance * 1000.0,
            (beam_x * pixel_x * 1000.0, beam_y * pixel_y * 1000.0),
            "+x",
            "-y",
            (1000 * pixel_x, 1000 * pixel_y),
            (nx, ny),
            (minimum_trusted_value, overload),
            [],
            px_mm=px_mm,
        )

        if material is not None:
            detector[0].set_thickness(thickness)
            detector[0].set_material(material)
            detector[0].set_mu(mu)

        return detector

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header, though for miniCBF
        there are limited options for this."""

        #  if "Phi" in self._cif_header_dictionary:
        #      phi_value = float(self._cif_header_dictionary["Phi"].split()[0])

        return self._goniometer_factory.single_axis()

    def _beam(self):
        """Return a simple model for the beam."""

        wavelength = float(self._cif_header_dictionary["Wavelength"].split()[0])

        beam = self._beam_factory.simple(wavelength)

        try:
            flux = float(self._cif_header_dictionary["Flux"].split()[0])
            beam.set_flux(flux)
        except (IndexError, KeyError):
            pass

        try:
            transmission = float(self._cif_header_dictionary["Transmission"].split()[0])
            beam.set_transmission(transmission)
        except (IndexError, KeyError):
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

        if cbf_header["byte_offset"]:
            pixel_values = uncompress(
                packed=data[data_offset : data_offset + cbf_header["size"]],
                fast=cbf_header["fast"],
                slow=cbf_header["slow"],
            )
        elif cbf_header["no_compression"]:
            assert len(self.get_detector()) == 1
            with self.open_file(self._image_file) as f:
                f.read(data_offset)
                pixel_values = read_int32(streambuf(f), cbf_header["length"])
            pixel_values.reshape(flex.grid(cbf_header["slow"], cbf_header["fast"]))

        else:
            raise ValueError(
                "Compression of type other than byte_offset or none is not supported (contact authors)"
            )

        return pixel_values

    def get_raw_data(self):
        if self._raw_data is None:
            data = self._read_cbf_image()
            self._raw_data = data

        return self._raw_data

    def detectorbase_start(self):
        self.detectorbase = PilatusImage(self._image_file)
        self.detectorbase.readHeader()  # necessary for LABELIT

    @staticmethod
    def as_file(
        detector,
        beam,
        gonio,
        scan,
        data,
        path,
        header_convention="GENERIC_MINI",
        det_type="GENERIC",
    ):
        """Note to developers: first attempt to write a miniCBF given a dxtbx-style experiment,
        But fields are not filled rigorously as in cbflib/src/cbf_minicbf_header.c
        Present code does not account for:
          Pilatus model number and serial number
          Data collection date and time
          Sensor material
          Dead time (exposure time - exposure period)
          Highest trusted value
          Energy threshold for photon counting
          Gain
          Bad pixels
          Auxiliary files (bad pixel mask, flat field, trim, image path)
          Detector not normal to beam
        """
        path = pathlib.Path(path)
        cbf = cbf_wrapper()
        cbf_root = path.with_suffix(".cbf")
        cbf.new_datablock(path.stem.encode())

        """Data items in the ARRAY_DATA category are the containers for the array data
        items described in the category ARRAY_STRUCTURE."""
        cbf.add_category("array_data", ["header_convention", "header_contents", "data"])

        # get pixel info out of detector object
        panel = detector[0]
        pixel_xy = panel.get_pixel_size()
        pixel_x_microns = pixel_xy[0] * 1000
        pixel_y_microns = pixel_xy[1] * 1000

        # make sure we get the right units
        thickness_meters = panel.get_thickness() / 1000.0

        material = panel.get_material()
        if material == "Si":
            sensor = "Silicon"
        elif material == "CdTe":
            sensor = "CdTe"
        elif panel.get_type() == "SENSOR_CCD":
            sensor = "CCD"
        else:
            sensor = "unknown"

        # maybe someday more people will do this
        flux = beam.get_flux()
        transmission = beam.get_transmission()
        polarization_fraction = beam.get_polarization_fraction()

        # get exposure information
        if scan is None:
            exposure_period = phi_start = osc = 0
        else:
            exposure_period = scan.get_exposure_times()[0]
            phi_start = scan.get_oscillation()[0]
            osc = scan.get_oscillation()[1]

        # automatically account for read-out time?
        # exposure_time = exposure_period-0.00203  # this is appropriate for Pilatus3 S
        exposure_time = exposure_period  # simulation is a perfect detector

        tau = 0  # assume simulation is a perfect detector with no pile-up error
        trusted_range = detector[0].get_trusted_range()
        count_cutoff = trusted_range[1]

        wavelength = beam.get_wavelength()  # get the wavelength in the conventional way
        energy = factor_ev_angstrom / wavelength
        threshold = energy / 2  # presume normal data collection convention

        bad_pixels = 0  # maybe get this from negative pixel values?

        assert len(detector) == 1, "only single-panel detectors supported"
        distance_meters = detector[0].get_distance() / 1000

        # fixed now, no longer backwards
        beam_center = detector[0].get_beam_centre_px(beam.get_s0())
        ORGX, ORGY = beam_center

        cbf.add_row(
            [
                header_convention,
                """
# Detector: %(det_type)s, S/N 60-0000
# 1972-01-01T00:00:00.000
# Pixel_size %(pixel_x_microns)ge-6 m x %(pixel_x_microns)ge-6 m
# %(sensor)s sensor, thickness %(thickness_meters)f m
# Exposure_time %(exposure_time).7f s
# Exposure_period %(exposure_period).7f s
# Tau = %(tau)f s
# Count_cutoff %(count_cutoff)d counts
# Threshold_setting: %(threshold)d eV
# Gain_setting: autog (vrf = 1.000)
# N_excluded_pixels = %(bad_pixels)d
# Excluded_pixels: badpix_mask.tif
# Flat_field: (nil)
# Trim_file: (nil)
# Image_path: /ramdisk/
# Wavelength %(wavelength)g A
# Detector_distance %(distance_meters)g m
# Beam_xy (%(ORGX)g, %(ORGY)g) pixels
# Flux %(flux)g
# Filter_transmission %(transmission).4f
# Start_angle %(phi_start).4f deg.
# Angle_increment %(osc).4f deg.
# Polarization %(polarization_fraction).3f
# Detector_2theta 0.0000 deg.
# Kappa 0.0000 deg.
# Phi %(phi_start).4f deg.
"""
                % locals(),
            ]
        )

        binary_id = 1
        focus = data.focus()
        data2 = data.copy_to_byte_str()
        elements = len(data)
        byteorder = b"little_endian"
        dimfast = focus[1]
        dimmid = focus[0]
        dimslow = 1
        padding = 0
        elsize = 4
        elsigned = 1

        cbf.set_integerarray_wdims_fs(
            pycbf.CBF_BYTE_OFFSET,
            binary_id,
            data2,
            elsize,
            elsigned,
            elements,
            byteorder,
            dimfast,
            dimmid,
            dimslow,
            padding,
        )

        cbf.write_widefile(
            str(cbf_root).encode(),
            pycbf.CBF,
            pycbf.MIME_HEADERS | pycbf.MSG_DIGEST | pycbf.PAD_4K,
            0,
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMini.understand(arg))
