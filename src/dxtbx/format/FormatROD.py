"""Support for the Rigaku Oxford Diffraction image format

NB: Rigaku datasets may use a non-zero-padded image incremental serial number.
At present, this is not compatible with assumptions in dxtbx. In order to
import these datasets, the images should be renumbered first.

See https://github.com/cctbx/dxtbx/issues/646 for details."""

# Takanori Nakane took David Waterman's code to parse headers from
#  https://github.com/cctbx/cctbx_project/commit/b95467f3b2a70a37eeb820ea294128a32551700c
# and heavily modified it. The original commit in the cctbx_project repository is orphan now.
from __future__ import annotations

__author__ = "David Waterman, Takanori Nakane"
__copyright__ = "Copyright 2018-2023 United Kingdom Research and Innovation & 2022-2023 Takanori Nakane"
__license__ = "BSD 3-clause"

import re
import struct

import numpy as np

from scitbx import matrix
from scitbx.array_family import flex
from scitbx.math import r3_rotation_axis_and_angle_as_matrix

from dxtbx.ext import uncompress_rod_TY6
from dxtbx.format.Format import Format
from dxtbx.model.beam import Probe
from dxtbx.model.detector import Detector


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
            try:
                hdr = f.read(256).decode("ascii")
            except UnicodeDecodeError:
                return False

        lines = hdr.splitlines()
        if len(lines) < 2:
            return False

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
            alpha2_wavelength = struct.unpack("<d", f.read(8))[0]
            alpha12_wavelength = struct.unpack("<d", f.read(8))[0]
            f.seek(offset + general_nbytes + special_nbytes + 640)
            # detector rotation in degrees along e1, e2, e3
            detector_rotns = struct.unpack("<ddd", f.read(24))
            # direct beam position when all angles are zero (FIXME: not completely sure)
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
            "alpha2_wavelength": alpha2_wavelength,
            "alpha12_wavelength": alpha12_wavelength,
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

        self._gonio_start_angles = (
            np.array(self._bin_header["start_angles_steps"])
            * np.nan_to_num(np.array(self._bin_header["step_to_rad"]))
            / np.pi
            * 180
        )
        self._gonio_end_angles = (
            np.array(self._bin_header["end_angles_steps"])
            * np.nan_to_num(np.array(self._bin_header["step_to_rad"]))
            / np.pi
            * 180
        )

        self._scan_axis = -1
        for axis in [0, 3]:  # 0 - OMEGA, 3 - PHI; the default is omega scan
            if self._gonio_start_angles[axis] != self._gonio_end_angles[axis]:
                self._scan_axis = axis
                break

        return

    # Rigaku/Oxford coordinate system:
    # - e3: parallel to OMEGA (and THETA), upwards
    # - e1: along the beam, from the crystal to the source
    # - e2: completes the right handed coordinate system
    #
    # The detector origin is at the lower left corner looking from the crystal (or the source).
    # - slow axis: (roughly) vertical and towards the ceiling
    # - fast axis: (roughly) horizontal and towards right (at THETA = 0)
    #
    # Importantly, positive THETA, OMEGA, PHI rotations are CLOCKWISE looking down from the ceiling,
    # while positive CHI=KAPPA rotation is CLOCKWISE looking from the detector at OMEGA=0, THETA=0.
    # See specifications uploaded in https://github.com/cctbx/dxtbx/issues/11#issue-434741559.
    #
    # In this dxtbx coordinate system:
    # - Z: along the beam, from the crystal to the source
    # - Y: vertical, towards the ceiling
    # - X: completes the right handed coordinate system
    # Thus, Z = e1, Y = e3, X = e2.
    #
    # The coordinate system of XDS.INP generated by CrysAlisPro:
    # - Z: along the beam, from the source to the crystal
    # - Y: vertical, towards the ground
    # - X: completes the right handed coordinate system
    # Thus, Y and Z are pointing the opposite directions from DIALS.

    def get_goniometer(self, index=None):
        return Format.get_goniometer(self)

    def _goniometer(self):
        if self._scan_axis == -1:
            raise NotImplementedError("Still shots not implemented yet.")
        elif self._scan_axis == 0:  # OMEGA
            dxtbx_scan_axis = 2
        elif self._scan_axis == 3:  # PHI
            dxtbx_scan_axis = 0
        else:
            pass  # should not happen

        # FIXME: sometimes XDS.INP generated by CrysAlisPro has
        #  tiny deviations from (0, 1, 0). I don't know how to calculate it.
        alpha_rad = self._bin_header["angles_in_deg"][0] * np.pi / 180
        axes = flex.vec3_double(
            ((0, -1, 0), (0, -np.cos(alpha_rad), np.sin(alpha_rad)), (0, -1, 0))
        )

        # angles[self._scan_axis] is not used anyway
        # angles are in degrees!
        angles = flex.double(
            (
                self._gonio_start_angles[3],
                self._gonio_start_angles[2],
                self._gonio_start_angles[0],
            )
        )
        names = flex.std_string(("PHI", "KAPPA=CHI", "OMEGA"))

        return self._goniometer_factory.make_multi_axis_goniometer(
            axes, angles, names, dxtbx_scan_axis
        )

    def get_beam(self, index=None):
        return Format.get_beam(self)

    def _beam(self):
        wavelength = self._bin_header["alpha12_wavelength"]
        if wavelength <= 0.05:
            probe = Probe.electron
        else:
            probe = Probe.xray

        return self._beam_factory.make_polarized_beam(
            sample_to_source=(0.0, 0.0, 1.0),
            wavelength=wavelength,
            polarization=(0, 1, 0),
            polarization_fraction=0.5,
            probe=probe,
        )

    def get_detector(self, index=None):
        return Format.get_detector(self)

    def _detector(self):
        theta_rad = self._gonio_start_angles[1] / 180 * np.pi
        detector_rotns_rad = np.array(self._bin_header["detector_rotns"]) / 180 * np.pi

        # I don't know why only rot_e1 is clockwise but this matches
        # DIRECTION_OF_DETECTOR_X-AXIS/Y-AXIS in XDS.INP from CrysAlisPro.
        # Note that XDS.INP's directions of Y and Z are opposite from ours.
        rot_e1 = np.array(
            r3_rotation_axis_and_angle_as_matrix([0, 0, 1], detector_rotns_rad[0])
        ).reshape(3, 3)  # clockwise along e1 = Z
        rot_e2 = np.array(
            r3_rotation_axis_and_angle_as_matrix([-1, 0, 0], detector_rotns_rad[1])
        ).reshape(3, 3)  # ANTI-clockwise along e2 = X
        rot_theta = np.array(
            r3_rotation_axis_and_angle_as_matrix([0, -1, 0], theta_rad)
        ).reshape(3, 3)  # ANTI-clockwise along e3 = Y
        detector_axes = rot_theta.dot(rot_e2.dot(rot_e1))

        pixel_size_x = self._bin_header["real_px_size_x"]
        pixel_size_y = self._bin_header["real_px_size_y"]
        origin_at_zero = np.array(
            [
                -self._bin_header["origin_px_x"] * pixel_size_x,
                -self._bin_header["origin_px_y"] * pixel_size_y,
                -self._bin_header["distance_mm"],
            ]
        )
        origin = detector_axes.dot(origin_at_zero)

        detector = self._detector_factory.make_detector(
            "PAD",
            detector_axes[:, 0],
            detector_axes[:, 1],
            origin,
            (pixel_size_x, pixel_size_y),
            (self._txt_header["NX"], self._txt_header["NY"]),
            (0, self._bin_header["overflow_threshold"]),
        )

        return detector

    def get_scan(self, index=None):
        return Format.get_scan(self)

    def _scan(self):
        """Return the scan information for this image."""

        start_angle = self._gonio_start_angles[self._scan_axis]
        end_angle = self._gonio_end_angles[self._scan_axis]

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
            return self._get_raw_data_ty6_native()
        else:
            raise NotImplementedError(f"Can't handle compression: {comp}")

    def _get_raw_data_ty6_native(self):
        offset = self._txt_header["NHEADER"]
        nx = self._txt_header["NX"]
        ny = self._txt_header["NY"]
        with open(self._image_file, "rb") as f:
            f.seek(offset)
            lbytesincompressedfield = struct.unpack("<l", f.read(4))[0]
            linedata = f.read(lbytesincompressedfield)
            offsets = f.read(4 * ny)

            return uncompress_rod_TY6(linedata, offsets, ny, nx)

    # Python implementation
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


class FormatROD_Arc(FormatROD):
    @staticmethod
    def understand(image_file):
        offset = 256
        general_nbytes = 512
        with FormatROD_Arc.open_file(image_file, "rb") as f:
            f.seek(offset + general_nbytes + 548)
            detector_type = struct.unpack("<l", f.read(4))[0]
            if detector_type in [12, 14]:
                return True
        return False

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
        bh = FormatROD._read_binary_header(
            image_file,
            offset,
            nbytes,
            general_nbytes,
            special_nbytes,
            km4gonio_nbytes,
            statistics_nbytes,
            history_nbytes,
        )

        # Seek to the end of the standard version 3 header and into the
        # extra camera parameters section
        with FormatROD_Arc.open_file(image_file, "rb") as f:
            f.seek(nbytes + 268)
            ix, iy, nx, ny, gapx, gapy = struct.unpack("<hhhhhh", f.read(12))
        assert ix in (2, 3)
        assert iy == 1
        assert nx == 385
        assert ny == 775
        assert gapx == 30
        assert gapy == 0
        bh["nx"] = nx
        bh["ny"] = ny
        bh["gap_px"] = gapx

        return bh

    def _detector(self):
        """2 or 3 panel detector, each rotated 38° from its neighbour."""

        theta_rad = self._gonio_start_angles[1] / 180 * np.pi
        detector_rotns_rad = np.array(self._bin_header["detector_rotns"]) / 180 * np.pi
        rot_e1 = np.array(
            r3_rotation_axis_and_angle_as_matrix([0, 0, 1], detector_rotns_rad[0])
        ).reshape(3, 3)  # clockwise along e1 = Z
        rot_e2 = np.array(
            r3_rotation_axis_and_angle_as_matrix([-1, 0, 0], detector_rotns_rad[1])
        ).reshape(3, 3)  # ANTI-clockwise along e2 = X
        rot_theta = np.array(
            r3_rotation_axis_and_angle_as_matrix([0, -1, 0], theta_rad)
        ).reshape(3, 3)  # ANTI-clockwise along e3 = Y
        detector_axes = rot_theta.dot(rot_e2.dot(rot_e1))

        pixel_size_x = self._bin_header["real_px_size_x"]
        pixel_size_y = self._bin_header["real_px_size_y"]
        origin_at_zero = np.array(
            [
                -self._bin_header["origin_px_x"] * pixel_size_x,
                -self._bin_header["origin_px_y"] * pixel_size_y,
                -self._bin_header["distance_mm"],
            ]
        )
        # Get origin for the first panel, not rotated
        origin = detector_axes.dot(origin_at_zero)

        d = Detector()

        root = d.hierarchy()
        root.set_local_frame(detector_axes[:, 0], detector_axes[:, 1], origin)

        self.coords = {}
        panel_idx = 0

        pnl_data = []
        nx = self._bin_header["nx"]
        ny = self._bin_header["ny"]
        gap_px = self._bin_header["gap_px"]
        pnl_data.append({"xmin": 0, "ymin": 0, "xmax": nx, "ymax": ny, "xmin_mm": 0})
        pnl_data.append(
            {
                "xmin": (nx + gap_px),
                "ymin": 0,
                "xmax": 2 * nx + gap_px,
                "ymax": ny,
                "xmin_mm": (nx + gap_px) * pixel_size_x,
            }
        )
        if self._bin_header["detector_type"] == 12:
            pnl_data.append(
                {
                    "xmin": (nx + gap_px) * 2,
                    "ymin": 0,
                    "xmax": 3 * nx + 2 * gap_px,
                    "ymax": ny,
                    "xmin_mm": (nx + gap_px) * 2 * pixel_size_x,
                }
            )

        # redefine fast, slow for the local frame
        fast = matrix.col((1.0, 0.0, 0.0))
        slow = matrix.col((0.0, 1.0, 0.0))

        for ipanel, pd in enumerate(pnl_data):
            xmin = pd["xmin"]
            xmax = pd["xmax"]
            ymin = pd["ymin"]
            ymax = pd["ymax"]
            xmin_mm = pd["xmin_mm"]

            origin_panel = fast * xmin_mm

            panel_name = "Panel%d" % panel_idx
            panel_idx += 1

            p = d.add_panel()
            p.set_type("SENSOR_PAD")
            p.set_name(panel_name)
            p.set_raw_image_offset((xmin, ymin))
            p.set_image_size((xmax - xmin, ymax - ymin))
            p.set_trusted_range((0, self._bin_header["overflow_threshold"]))
            p.set_pixel_size((pixel_size_x, pixel_size_y))
            p.set_local_frame(fast.elems, slow.elems, origin_panel.elems)
            p.set_projection_2d((-1, 0, 0, -1), (0, pd["xmin"]))
            self.coords[panel_name] = (xmin, ymin, xmax, ymax)

        # Now rotate the external panels. For the time being do the rotation around
        # and axis along the centre of the gap between the panels.
        if len(pnl_data) == 2:
            angle = 19.0
        elif len(pnl_data) == 3:
            angle = 38.0
        else:
            raise ValueError("Unexpected number of panels")

        left_panel = d[0]
        right_panel = d[len(d) - 1]

        # Point to rotate the left panel around, from the local origin
        pnl_fast = matrix.col(left_panel.get_local_fast_axis())
        pnl_slow = matrix.col(left_panel.get_local_slow_axis())
        pt = pnl_fast * (
            left_panel.get_image_size_mm()[0] + (gap_px / 2) * pixel_size_x
        )

        rotated = (-1.0 * pt).rotate_around_origin(pnl_slow, angle, deg=True)
        new_origin = pt + rotated
        new_fast = pnl_fast.rotate_around_origin(pnl_slow, angle, deg=True)
        left_panel.set_local_frame(new_fast.elems, pnl_slow.elems, new_origin.elems)

        # Point to rotate the right panel around, from the panel's origin
        pnl_fast = matrix.col(right_panel.get_local_fast_axis())
        pnl_slow = matrix.col(right_panel.get_local_slow_axis())
        pnl_origin = matrix.col(right_panel.get_local_origin())
        pt = -pnl_fast * (gap_px / 2) * pixel_size_x

        rotated = (-1.0 * pt).rotate_around_origin(pnl_slow, -angle, deg=True)
        new_origin = pnl_origin + pt + rotated
        new_fast = pnl_fast.rotate_around_origin(pnl_slow, -angle, deg=True)
        right_panel.set_local_frame(new_fast.elems, pnl_slow.elems, new_origin.elems)

        return d

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        raw_data = super().get_raw_data()

        # split into separate panels
        self._raw_data = []
        d = self.get_detector()
        for panel in d:
            xmin, ymin, xmax, ymax = self.coords[panel.get_name()]
            self._raw_data.append(raw_data[ymin:ymax, xmin:xmax])

        return tuple(self._raw_data)


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        if FormatROD.understand(arg):
            reader = FormatROD(arg)
            print(reader._txt_header)
            print(reader._bin_header)
            print()
            print(
                "Starting angles in degrees (OMEGA, THETA, KAPPA=CHI, PHI, OMEGA PRIME, THETA PRIME)"
            )
            print(reader._gonio_start_angles)
            print("Ending angles in degrees")
            print(reader._gonio_end_angles)
        else:
            print("Unsupported format")
