from __future__ import annotations

import math
import os
import re
import sys
import time

from boost_adaptbx.boost.python import streambuf
from scitbx.array_family import flex

from dxtbx import IncorrectFormatError
from dxtbx.ext import is_big_endian, read_uint16, read_uint16_bs
from dxtbx.format.Format import Format


class FormatNoniusKappaCCD(Format):
    """A class for reading files produced by 1990s era Nonius Kappa CCD instruments."""

    @staticmethod
    def read_header_lines(image_path):
        """
        Adapted from Fabio KcdImage _readheader
        """

        hdr_lines = []
        end_of_headers = False
        linect = 0
        with Format.open_file(image_path, "rb") as f:
            while not end_of_headers:
                one_line = f.readline()
                linect += 1
                try:
                    one_line = one_line.decode("ASCII")
                except UnicodeDecodeError:
                    end_of_headers = True
                else:
                    if len(one_line) > 100:
                        end_of_headers = True
                if not end_of_headers:
                    if one_line.strip() == "Binned mode":
                        one_line = "Mode = Binned"
                    if "=" not in one_line and linect > 1:
                        end_of_headers = True
                    if not end_of_headers:
                        hdr_lines.append(one_line)

        return hdr_lines

    @staticmethod
    def parse_header(header_lines):
        header_dic = {}

        for l in header_lines:
            separator = l.find("=")
            if separator == -1:
                continue
            k = l[:separator].strip()
            v = l[separator + 1 :].strip()
            if k in header_dic:
                header_dic[k] = header_dic[k] + "\n" + v
            else:
                header_dic[k] = v

        return header_dic

    @staticmethod
    def understand(image_file):
        """
        Look for characteristic character sequences in the
        header.
        """
        try:
            with Format.open_file(image_file, "rb") as fh:
                tag = fh.read(1024).decode("latin-1", "replace")
        except OSError:
            return False
        if tag[0:2] != "No":
            return False
        if re.search("^Kappa-support angle", tag, re.MULTILINE) is None:
            return False
        if re.search("^Binned mode", tag, re.MULTILINE) is None:
            return False
        if re.search("^Shutter = Closed", tag, re.MULTILINE) is not None:
            return False  # ignore dark current measurement
        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        super().__init__(str(image_file), **kwargs)

    def _start(self):
        self.headers = self.parse_header(self.read_header_lines(self._image_file))

    def _goniometer(self):
        """Return model for a kappa goniometer"""
        alpha = float(self.headers["Kappa-support angle"])
        omega = float(self.headers["Omega start"])
        kappa = float(self.headers["Kappa start"])
        phi = float(self.headers["Phi start"])
        direction = "-z"  #

        if "Phi scan range" in self.headers:
            scan_axis = "phi"
        else:
            scan_axis = "omega"

        return self._goniometer_factory.make_kappa_goniometer(
            alpha, omega, kappa, phi, direction, scan_axis
        )

    def _detector(self):
        """It appears that pixel size reported in the header does not
        account for binning, which doubles the size in both dimensions."""

        pix_fac = 0.001
        if self.headers["Mode"] == "Binned":
            pix_fac = 0.002

        pixel_size = [
            float(self.headers["pixel X-size (um)"]) * pix_fac,
            float(self.headers["pixel Y-size (um)"]) * pix_fac,
        ]

        image_size = [
            float(self.headers["X dimension"]),
            float(self.headers["Y dimension"]),
        ]

        beam_centre = (
            image_size[0] * pixel_size[0] * 0.5,
            image_size[1] * pixel_size[1] * 0.5,
        )

        gain = self.headers.get("Detector gain (ADU/count)", "1.0")
        gain = float(gain)

        return self._detector_factory.two_theta(
            sensor="CCD",
            distance=float(self.headers["Dx start"]),
            beam_centre=beam_centre,
            fast_direction="+y",  # Increasing to the right
            slow_direction="+x",  # Down, same as rotation axis
            two_theta_direction="+x",  # Same sense as omega and phi
            two_theta_angle=float(self.headers["Theta start"]) * 2.0,
            pixel_size=pixel_size,
            image_size=image_size,
            gain=gain,
            identifier=self.headers["Detector ID"],
            trusted_range=(0.0, 131070.0),  # May readout UInt16 twice
        )

    def _beam(self):
        """Return a simple model for a lab-based beam. As dxtbx
        has no laboratory target model, take the weighted values
        of the alpha1/alpha2 wavelengths."""

        a1 = float(self.headers["Alpha1"])
        a2 = float(self.headers["Alpha2"])
        wt = float(self.headers["Alpha1/Alpha2 ratio"])
        wavelength = (wt * a1 + a2) / (wt + 1.0)
        direction = [0.0, 0.0, 1.0]
        polarisation, pol_dir = self.get_beam_polarisation()
        return self._beam_factory.complex(
            sample_to_source=direction,
            wavelength=wavelength,
            polarization_plane_normal=pol_dir,
            polarization_fraction=polarisation,
        )

    def _scan(self):
        """All scan dates will be in the range 1969-2068 as per
        Python strptime. Please retire your CCD before 2069."""

        if "Phi scan range" in self.headers:
            rotax = "Phi"
        else:
            rotax = "Omega"

        osc_start = float(self.headers["%s start" % rotax])
        osc_range = float(self.headers["%s scan range" % rotax])
        exposure_time = float(self.headers["Exposure time"])

        simpletime = self.headers["Date"]
        epoch = time.mktime(time.strptime(simpletime, "%a %m/%d/%y    %I:%M:%S %p"))

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )

    def get_beam_polarisation(self):
        """Polarisation for single reflection off graphite 002
        monochromator as per manual. Hard-coded angles for Mo, Cu and
        Ag targets."""

        if self.headers["Target material"] == "MO":
            two_theta = 12.2  # From manual
        elif self.headers["Target material"] == "CU":
            two_theta = 26.6  # From manual
        elif self.headers["Target material"] == "AG":
            two_theta = 9.62  # Calculated

        # Assume an ideally imperfect monochromator, with the axis
        # of rotation of the monochromator in the horizontal direction
        # in a laboratory setup.

        pol_frac = 1.0 / (1 + math.cos(two_theta * math.pi / 180) ** 2)
        pol_dir = [0.0, 1.0, 0.0]
        return pol_frac, pol_dir

    def get_raw_data(self):
        """Return raw data from a Nonius CCD frame."""

        # Frame file contains a series of readouts, each pixel is
        # an unsigned little-endian 16-bit integer

        dim1 = int(self.headers["X dimension"])
        dim2 = int(self.headers["Y dimension"])
        nbReadOut = int(self.headers["Number of readouts"])
        expected_size = dim1 * dim2 * 2 * nbReadOut

        # Not clear if the image offset is present in the header,
        # therefore we seek from the end of the file.

        with self.open_file(self._image_file, "rb") as infile:
            fileSize = os.stat(self._image_file)[6]
            infile.seek(fileSize - expected_size, os.SEEK_SET)

            for i in range(nbReadOut):
                if is_big_endian():
                    raw_data = read_uint16_bs(streambuf(infile), dim1 * dim2)
                else:
                    raw_data = read_uint16(streambuf(infile), dim1 * dim2)

                if i == 0:
                    data = raw_data
                else:
                    data += raw_data

        data.reshape(flex.grid(dim2, dim1))
        return data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatNoniusKappaCCD.understand(arg))
