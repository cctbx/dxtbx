"""Support for the Rigaku Oxford Diffraction image format"""

# Takanori Nakane took David Waterman's code to parse headers from
#  https://github.com/cctbx/cctbx_project/commit/b95467f3b2a70a37eeb820ea294128a32551700c
# and heavily modified it. The original commit in the cctbx_project repository is orphan now.

from __future__ import annotations

__author__ = "David Waterman, Takanori Nakane"
__copyright__ = (
    "Copyright 2018 United Kingdom Research and Innovation & 2022 Takanori Nakane"
)
__license__ = "BSD 3-clause"

import re
import struct

import numpy as np

from scitbx.array_family import flex

from dxtbx.format.Format import Format


class FormatROD(Format):
    """An image reading class for Rigaku Oxford Diffraction images.

    The Rigaku Oxford Diffraction format is used by CrysAlis(Pro).
    """

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        Format.__init__(self, image_file, **kwargs)

    @staticmethod
    def understand(image_file):

        with FormatROD.open_file(image_file, "rb") as f:
            hdr = f.read(256).decode("ascii")
        lines = hdr.splitlines()

        vers = lines[0].split()
        if len(vers) < 2 or vers[0] != "OD" or vers[1] != "SAPPHIRE":
            return False

        compression = lines[1].split("=")
        if compression[0] != "COMPRESSION":
            return False

        return True

    @staticmethod
    def _read_ascii_header(image_file):
        """Read the ASCII header comprising the first 256 bytes of the file"""

        hd = {}
        with FormatROD.open_file(image_file, "rb") as f:
            hdr = f.read(256).decode("ascii")
        lines = hdr.splitlines()

        vers = lines[0].split()
        if len(vers) < 2 or vers[0] != "OD" or vers[1] != "SAPPHIRE":
            raise ValueError("Wrong header format")
        hd["version"] = float(vers[-1])

        compression = lines[1].split("=")
        if compression[0] != "COMPRESSION":
            raise ValueError("Wrong header format")
        hd["compression"] = compression[1]

        # Extract definitions from the 3rd - 5th line
        defn = re.compile(r"([A-Z]+=[ 0-9]+)")
        for line in lines[2:5]:
            sizes = defn.findall(line)
            for s in sizes:
                n, v = s.split("=")
                hd[n] = int(v)

        hd["time"] = lines[5].split("TIME=")[-1].strip("\x1a").rstrip()

        return hd

    @staticmethod
    def _read_binary_header(
        image_file,
        offset=256,
        nbytes=5120,
        general_nbytes=512,
        special_nbytes=768,
        km4gonio_nbytes=1024,
        statistics_nbytes=512,
        history_nbytes=2048,
    ):
        """Read the most relevant parameters from the binary header"""

        with FormatROD.open_file(image_file, "rb") as f:
            # General section
            f.seek(offset)
            bin_x, bin_y = struct.unpack("<hh", f.read(4))
            f.seek(offset + 22)
            chip_npx_x, chip_npx_y, im_npx_x, im_npx_y = struct.unpack(
                "<hhhh", f.read(8)
            )
            f.seek(offset + 36)
            num_points = struct.unpack("<I", f.read(4))[0]
            if num_points != im_npx_x * im_npx_y:
                raise ValueError("Cannot interpret binary header")

            # Special section
            f.seek(offset + general_nbytes + 56)
            gain = struct.unpack("<d", f.read(8))[0]
            f.seek(offset + general_nbytes + 464)
            overflow_flag, overflow_after_remeasure_flag = struct.unpack(
                "<hh", f.read(4)
            )
            f.seek(offset + general_nbytes + 472)
            overflow_threshold = struct.unpack("<l", f.read(4))[0]
            f.seek(offset + general_nbytes + 480)
            exposure_time_sec, overflow_time_sec = struct.unpack("<dd", f.read(16))
            f.seek(offset + general_nbytes + 548)
            detector_type = struct.unpack("<l", f.read(4))[0]
            f.seek(offset + general_nbytes + 568)
            real_px_size_x, real_px_size_y = struct.unpack("<dd", f.read(16))

            # Goniometer section
            f.seek(offset + general_nbytes + special_nbytes + 284)
            # angles for OMEGA, THETA, CHI(=KAPPA), PHI,
            # OMEGA_PRIME (also called DETECTOR_AXIS; what's this?), THETA_PRIME
            start_angles_steps = struct.unpack("<llllllllll", f.read(40))
            end_angles_steps = struct.unpack("<llllllllll", f.read(40))
            f.seek(offset + general_nbytes + special_nbytes + 368)
            step_to_rad = struct.unpack("<dddddddddd", f.read(80))
            f.seek(offset + general_nbytes + special_nbytes + 552)
            # FIXME: I don't know what these are. Isn't the beam along e1 by definition??
            beam_rotn_around_e2, beam_rotn_around_e3 = struct.unpack("<dd", f.read(16))
            alpha1_wavelength = struct.unpack("<d", f.read(8))[0]
            f.seek(offset + general_nbytes + special_nbytes + 640)
            # detector rotation in degrees along e1, e2, e3
            detector_rotns = struct.unpack("<ddd", f.read(24))
            # FIXME: direct beam position when all angles are zero?
            origin_px_x, origin_px_y = struct.unpack("<dd", f.read(16))
            # alpha and beta are angles between KAPPA(=CHI) and THETA, and e3.
            angles_in_deg = struct.unpack(
                "<dddd", f.read(32)
            )  # alpha, beta, gamma, delta
            f.seek(offset + general_nbytes + special_nbytes + 712)
            distance_mm = struct.unpack("<d", f.read(8))[0]

        return {
            "bin_x": bin_x,
            "bin_y": bin_y,
            "chip_npx_x": chip_npx_x,
            "chip_npx_y": im_npx_y,
            "im_npx_x": im_npx_x,
            "im_npx_y": im_npx_y,
            "gain": gain,
            "overflow_flag": overflow_flag,
            "overflow_after_remeasure_flag": overflow_after_remeasure_flag,
            "overflow_threshold": overflow_threshold,
            "exposure_time_sec": exposure_time_sec,
            "overflow_time_sec": overflow_time_sec,
            "detector_type": detector_type,
            "real_px_size_x": real_px_size_x,
            "real_px_size_y": real_px_size_y,
            "start_angles_steps": start_angles_steps,
            "end_angles_steps": end_angles_steps,
            "step_to_rad": step_to_rad,
            "beam_rotn_around_e2": beam_rotn_around_e2,
            "beam_rotn_around_e3": beam_rotn_around_e3,
            "alpha1_wavelength": alpha1_wavelength,
            "detector_rotns": detector_rotns,
            "origin_px_x": origin_px_x,
            "origin_px_y": origin_px_y,
            "angles_in_deg": angles_in_deg,
            "distance_mm": distance_mm,
        }

    def _start(self):
        """Open the image file, read useful metadata into an internal dictionary
        self._header_dictionary"""

        self._txt_header = self._read_ascii_header(self._image_file)
        self._bin_header = self._read_binary_header(self._image_file)

        return

    # Rigaku/Oxford Geometry:
    # - e3: parallel to OMEGA (and THETA)
    # - e1: crystal to source
    # - e2: completes the right handed coordinate system

    def get_goniometer(self, index=None):
        return Format.get_goniometer(self)

    def _goniometer(self):
        # FIXME: sometimes XDS.INP generated by CrysAlisPro has
        #  tiny deviationis from (0, 1, 0). I don't know how to calculate it.
        # DIALS' third axis points to the opposite of XDS; thus the direction will be opposite.
        direction = (0.0, -1.0, 0.0)

        # FIXME: represent kappa axis in fixed_rotation
        return self._goniometer_factory.known_axis(direction)

    def get_beam(self, index=None):
        return Format.get_beam(self)

    def _beam(self):
        return self._beam_factory.make_polarized_beam(
            sample_to_source=(0.0, 0.0, 1.0),
            wavelength=self._bin_header["alpha1_wavelength"],
            polarization=(0, 1, 0),
            polarization_fraction=0.5,
        )

    def get_detector(self, index=None):
        return Format.get_detector(self)

    def _detector(self):
        gonio_angles = np.array(self._bin_header["start_angles_steps"]) * np.array(
            self._bin_header["step_to_rad"]
        )
        detector_rotns_rad = np.array(self._bin_header["detector_rotns"]) / 180 * np.pi
        rot_e1 = np.array(
            [
                np.cos(detector_rotns_rad[0]),
                np.sin(detector_rotns_rad[0]),
                0,
                -np.sin(detector_rotns_rad[0]),
                np.cos(detector_rotns_rad[0]),
                0,
                0,
                0,
                1,
            ]
        ).reshape(3, 3)
        rot_e2 = np.array(
            [
                1,
                0,
                0,
                0,
                np.cos(detector_rotns_rad[1]),
                np.sin(detector_rotns_rad[1]),
                0,
                -np.sin(detector_rotns_rad[1]),
                np.cos(detector_rotns_rad[1]),
            ]
        ).reshape(3, 3)
        rot_theta = np.array(
            [
                np.cos(gonio_angles[1]),
                0,
                np.sin(gonio_angles[1]),
                0,
                1,
                0,
                -np.sin(gonio_angles[1]),
                0,
                np.cos(gonio_angles[1]),
            ]
        ).reshape(3, 3)
        detector_axes = rot_theta.dot(rot_e2.dot(rot_e1))
        # The third axis points to the opposite of XDS.
        detector_axes[2, :] *= -1

        pixel_size_x = self._bin_header["real_px_size_x"]
        pixel_size_y = self._bin_header["real_px_size_y"]
        origin_at_zero = np.array(
            [
                -self._bin_header["origin_px_x"] * pixel_size_x,
                +self._bin_header["origin_px_y"] * pixel_size_y,
                +self._bin_header["distance_mm"],
            ]
        )

        # FIXME: this formula seem to give the right answer but I don't know why.
        # XDS exporter in IPR's CrysAlisPro seems broken. It writes the same ORGX/Y
        # regardless of the theta angles.
        origin = detector_axes.dot(origin_at_zero)

        detector = self._detector_factory.make_detector(
            "PAD",
            detector_axes[:, 0],
            -detector_axes[:, 1],
            origin,
            (pixel_size_x, pixel_size_y),
            (self._txt_header["NX"], self._txt_header["NY"]),
            (0, self._bin_header["overflow_threshold"]),
        )  # not sure about min

        return detector

    def get_scan(self, index=None):
        return Format.get_scan(self)

    def _scan(self):
        """Return the scan information for this image."""

        for axis in [0, 3]:  # 0 - OMEGA, 3 - PHI: FIXME: is it always PHI scan?
            start_angle = (
                self._bin_header["start_angles_steps"][axis]
                * self._bin_header["step_to_rad"][axis]
                / np.pi
                * 180
            )
            end_angle = (
                self._bin_header["end_angles_steps"][axis]
                * self._bin_header["step_to_rad"][axis]
                / np.pi
                * 180
            )
            if start_angle != end_angle:
                break

        return self._scan_factory.single_file(
            filename=self._image_file,
            exposure_times=self._bin_header["exposure_time_sec"],
            osc_start=start_angle,
            osc_width=end_angle - start_angle,
            epoch=None,
        )

    def get_raw_data(self):
        comp = self._txt_header["compression"].strip()
        if comp.startswith("TY6"):
            return self._get_raw_data_ty6()
        else:
            raise NotImplementedError("Can't handle compression: {0}".format(comp))

    def _get_raw_data_ty6(self):
        offset = self._txt_header["NHEADER"]
        nx = self._txt_header["NX"]
        ny = self._txt_header["NY"]
        with open(self._image_file, "rb") as f:
            f.seek(offset)
            lbytesincompressedfield = struct.unpack("<l", f.read(4))[0]
            linedata = np.fromfile(f, dtype=np.uint8, count=lbytesincompressedfield)
            offsets = struct.unpack("<%dI" % ny, f.read(4 * ny))

            image = np.zeros((ny, nx), dtype=np.int32)
            for iy in range(ny):
                image[iy, :] = self.decode_TY6_oneline(linedata[offsets[iy] :], nx)
            return flex.int(image)

    def decode_TY6_oneline(self, linedata, w):
        """Decompress TY6 encoded pixels for a single line.
        w is the number of pixels in the fast axis."""

        BLOCKSIZE = 8
        SHORT_OVERFLOW = np.int32(254)
        LONG_OVERFLOW = 255
        SHORT_OVERFLOW_SIGNED = SHORT_OVERFLOW - 127
        LONG_OVERFLOW_SIGNED = LONG_OVERFLOW - 127

        ipos = 0
        opos = 0
        ret = np.zeros(w, dtype=np.int32)

        nblock = (w - 1) // (BLOCKSIZE * 2)
        nrest = (w - 1) % (BLOCKSIZE * 2)

        firstpx = linedata[ipos]
        ipos += 1
        if firstpx < SHORT_OVERFLOW:
            ret[opos] = firstpx - 127
        elif firstpx == LONG_OVERFLOW:
            ret[opos] = linedata[ipos : (ipos + 4)].view(np.int32)[0]
            ipos += 4
        else:
            ret[opos] = linedata[ipos : (ipos + 2)].view(np.int16)[0]
            ipos += 2
        opos += 1

        for k in range(nblock):
            bittype = linedata[ipos]
            nbits = (bittype & 15, (bittype >> 4) & 15)
            ipos += 1
            # ipos_bit = ipos * 8

            for i in range(2):
                nbit = nbits[i]
                # 0: 0
                # 1: 1 or 0
                # 2: -1, 0, 1, 2
                # 3: -3, -2, -1, 0, 1, 2, 3, 4 (= b111 -3 = 7 - 3)
                zero_at = 0
                if nbit > 1:
                    zero_at = (1 << (nbit - 1)) - 1

                v = 0
                for j in range(nbit):
                    v |= linedata[ipos] << (8 * j)
                    ipos += 1

                mask = (1 << nbit) - 1
                for j in range(BLOCKSIZE):
                    ret[opos] = ((v >> (nbit * j)) & mask) - zero_at
                    opos += 1

            for i in range(opos - BLOCKSIZE * 2, opos):
                offset = ret[i]

                # surprisingly, this IF is very slow!
                if offset >= SHORT_OVERFLOW_SIGNED:
                    if offset >= LONG_OVERFLOW_SIGNED:
                        offset = linedata[ipos : (ipos + 4)].view(np.int32)[0]
                        ipos += 4
                    else:
                        offset = linedata[ipos : (ipos + 2)].view(np.int16)[0]
                        ipos += 2

                ret[i] = offset + ret[i - 1]

        for i in range(nrest):
            px = linedata[ipos]
            ipos += 1
            if px < SHORT_OVERFLOW:
                ret[opos] = ret[opos - 1] + px - 127
            elif px == LONG_OVERFLOW:
                ret[opos] = (
                    ret[opos - 1] + linedata[ipos : (ipos + 4)].view(np.int32)[0]
                )
                ipos += 4
            else:
                ret[opos] = (
                    ret[opos - 1] + linedata[ipos : (ipos + 2)].view(np.int16)[0]
                )
                ipos += 2
            opos += 1

        return ret


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        if FormatROD.understand(arg):
            reader = FormatROD(arg)
            print(reader._txt_header)
            print(reader._bin_header)
        else:
            print("Unsupported format")
