from __future__ import annotations

import io
import math
import os
import re
import sys
import time

from boost_adaptbx.boost.python import streambuf
from scitbx.array_family import flex

from dxtbx import IncorrectFormatError
from dxtbx.ext import read_uint16_bs
from dxtbx.format.Format import Format

MINIMUM_KEYS = [
    # 'ByteOrder', Assume little by default
    "Data type",
    "X dimension",
    "Y dimension",
    "Number of readouts",
]


class FormatNoniusKappaCCD(Format):
    """This format produced by 1990s era Nonius Kappa CCD instruments.
    Consists of a header containing lines of the form XXX = YYY, followed by an image.
    We reuse Fabio code.
    """

    @staticmethod
    def read_header_lines(image_path):
        """
        Adapted from Fabio KcdImage _readheader
        """

        hdr_lines = []
        end_of_headers = False
        with Format.open_file(image_path, "rb") as f:
            while not end_of_headers:
                one_line = f.readline()
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
                    if "=" not in one_line:
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
        Try to find some characteristic character sequences in the
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

        # Kappa axis points into the kappa head as clockwise rotation when
        # viewed from above. The datum position has the kappa head under
        # the incoming beam. All primary axes rotate clockwise when viewed
        # from above, so the kappa axis has positive x, z components,
        # corresponding to "+z" when the calculation in in make_kappa_goniometer
        # is reviewed.

        direction = "+z"  #

        if "Phi scan range" in self.headers:
            scan_axis = "phi"
        else:
            scan_axis = "omega"

        return self._goniometer_factory.make_kappa_goniometer(
            alpha, omega, kappa, phi, direction, scan_axis
        )

    def _detector(self):
        """It appears that pixel size reported in the header does not
        account for binning, which doubles the size. CHECK"""

        pix_fac = 0.001
        # if self.headers["Mode"] == "Binned":
        #    pix_fac = 0.002

        pixel_size = [
            float(self.headers["pixel X-size (um)"]) * pix_fac,
            float(self.headers["pixel Y-size (um)"]) * pix_fac,
        ]

        image_size = [
            float(self.headers["X dimension"]),
            float(self.headers["Y dimension"]),
        ]

        self._detector_factory.two_theta(
            sensor="CCD",
            distance=float(self.headers["Dx start"]),
            beam_centre=[image_size[0] / 2.0, image_size[1] / 2],
            fast_direction="+y",  # Increasing to the right
            slow_direction="+x",  # Down, same as rotation axis
            two_theta_direction=[1, 0, 0],  # Same sense as omega and phi
            two_theta_angle=float(self.headers["Theta start"]) * 2.0,
            pixel_size=pixel_size,
            image_size=image_size,
            identifier=self.headers["Detector ID"],
        )

    def _beam(self):
        """Return a simple model for a lab-based beam. As dxtbx
        has no laboratory target model, take the weighted values
        of the alpha1/alpha2 wavelengths. Polarisation is that
        obtained after single reflection from graphite 002
        monochromator."""
        a1 = float(self.headers["Alpha1"])
        a2 = float(self.headers["Alpha2"])
        wt = float(self.headers["Alpha1/Alpha2 ratio"])
        wavelength = wt * a1 + a2
        direction = [0.0, 0.0, 1.0]  # imgCIF standard
        polarisation, pol_dir = self.get_beam_polarisation()
        return self._beam_factory.complex(
            sample_to_source=direction,
            wavelength=wavelength,
            polarization=pol_dir,
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

        simpletime = self.header["Date"]
        epoch = time.mktime(time.strptime(simpletime, "%a %m/%d/%y    %I:%M:%S %p"))

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )

    def get_beam_polarisation(self):
        """Polarisation for single reflection off graphite
        002 monochromator. Hard-coded angles for MO and CU
        targets
        """
        if self.headers["Target material"] == "MO":
            two_theta = 12.2  # From manual
        elif self.headers["Target material"] == "CU":
            two_theta = 26.6  # From manual
        elif self.headers["Target material"] == "AG":
            two_theta = 9.62  # Calculated

        pol_frac = (1 + math.cos(two_theta * math.pi / 180) ** 2) / 2
        pol_dir = [1.0, 0.0, 0.0]  # To be confirmed later.
        return pol_frac, pol_dir

    def get_raw_data(self):
        # Compute image size: always little-endian UInt16 = 2 bytes

        dim1 = int(self.headers["X dimension"])
        dim2 = int(self.headers["Y dimension"])
        nbReadOut = int(self.headers["Number of readouts"])
        expected_size = dim1 * dim2 * 2 * nbReadOut

        # Now seek to beginning counting from the end

        with self.open_file(self._image_file, "rb") as infile:
            try:
                infile.seek(-expected_size, io.SEEK_END)
            except Exception:
                print(
                    "Warning: Seeking from end not implemented for file %s"
                    % self._image_file
                )
                if hasattr(infile, "measure_size"):
                    fileSize = infile.measure_size()
                elif hasattr(infile, "size"):
                    fileSize = infile.size
                elif hasattr(infile, "getSize"):
                    fileSize = infile.getSize()
                else:
                    print("Unable to guess the file-size of %s" % self._image_file)
                    fileSize = os.stat(self._image_file)[6]
                infile.seek(fileSize - expected_size - infile.tell(), 1)

        # Read data into array. Data are little endian based on Fabio approach
        # Data may be double-read, presumably in case some CCD bins had something
        # left after the first read.

        for i in range(nbReadOut):
            raw_data = read_uint16_bs(streambuf(infile), dim1 * dim2)
            if i == 0:
                data = raw_data
            else:
                data += raw_data

        data.reshape(flex.grid(dim2, dim1))  # Not sure about correct order of dims
        return data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatNoniusKappaCCD.understand(arg))
