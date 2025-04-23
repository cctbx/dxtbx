from __future__ import annotations

from boost_adaptbx.boost.python import streambuf
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx import IncorrectFormatError
from dxtbx.ext import (
    is_big_endian,
    read_uint8,
    read_uint16,
    read_uint16_bs,
    read_uint32,
    read_uint32_bs,
)
from dxtbx.format.FormatBruker import FormatBruker
from dxtbx.model import SimplePxMmStrategy
from dxtbx.model.beam import Probe


class FormatBrukerELA(FormatBruker):
    @staticmethod
    def understand(image_file):
        try:
            header_lines = FormatBruker.read_header_lines(image_file)
        except OSError:
            return False

        header_dic = FormatBruker.parse_header(header_lines)

        dettype = header_dic.get("DETTYPE").upper()
        if dettype is None:
            return False

        if "DECTRIS-ELA" not in dettype:
            return False

        return True

    def _start(self):
        try:
            header_lines = FormatBruker.read_header_lines(self._image_file)
        except OSError:
            return False

        self.header_dict = FormatBrukerELA.parse_header(header_lines)

        # The ELA format can't currently use BrukerImage, for the same reason
        # as for the Photon-II/III (hardcoded image size is one factor)
        # https://github.com/cctbx/cctbx_project/issues/65
        # from iotbx.detectors.bruker import BrukerImage
        # self.detectorbase = BrukerImage(self._image_file)

    def _goniometer(self):
        # goniometer angles in ANGLES are 2-theta, omega, phi, chi (FIXED)
        # AXIS indexes into this list to define the scan axis (in FORTRAN counting)
        # START and RANGE define the start and step size for each image

        _, omega, phi, chi = map(float, self.header_dict["ANGLES"].split())
        scan_axis = ["NONE", "2THETA", "OMEGA", "PHI", "CHI", "X", "Y", "Z"]
        scan_axis = scan_axis[int(self.header_dict["AXIS"])]
        names = flex.std_string(("PHI", "CHI", "OMEGA"))
        scan_axis = flex.first_index(names, scan_axis)
        if scan_axis is None:
            scan_axis = 2  # "OMEGA" default

        # Axes here may be incorrect
        axes = flex.vec3_double(((1, 0, 0), (0, 0, -1), (0, 1, 0)))
        omega -= 180
        angles = flex.double((phi, chi, omega))

        g = self._goniometer_factory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis
        )

        # It is preferable to return a single axis goniometer, as this can be
        # further optimised by dials.find_rotation_axis
        return self._goniometer_factory.known_axis(g.get_rotation_axis())

    def _calculate_gain(self, wavelength):
        """The CCDPARM header item contains 5 items. For an X-ray detector
        are described by:
            1. readnoise
            2. e/ADU
            3. e/photon
            4. bias
            5. full scale
        and the gain in ADU/X-ray is given by (e/photon) / (e/ADU).
        It is not exactly clear what the gain is for this electron detector, but
        these files seem to set it to 1
        """
        ccdparm = self.header_dict["CCDPARM"].split()
        e_ADU = float(ccdparm[1])
        e_photon = float(ccdparm[2])
        if e_ADU == 0:
            return 1.0
        gain = e_photon / e_ADU
        return gain

    def _detector(self):
        # goniometer angles in ANGLES are 2-theta, omega, phi, chi (FIXED)
        two_theta = float(self.header_dict["ANGLES"].split()[0])

        # Assume Bruker full_scale value means saturation
        full_scale = float(self.header_dict["CCDPARM"].split()[-1])
        min_trusted_value = 0

        fast = matrix.col((1, 0, 0))
        slow = matrix.col((0, -1, 0))
        beam = matrix.col((0, 0, 1))
        pixel_mm = 0.075
        beam_pixel = [float(bp) for bp in self.header_dict["CENTER"].split()[:-3:-1]]
        distance_mm = 10.0 * float(self.header_dict["DISTANC"].split()[1])
        origin = (
            -distance_mm * beam
            - fast * pixel_mm * beam_pixel[1]
            - slow * pixel_mm * beam_pixel[0]
        )
        # 2theta rotation appears to be around the slow axis
        origin = origin.rotate_around_origin(slow, two_theta, deg=True)
        fast = fast.rotate_around_origin(slow, two_theta, deg=True)
        pixel_size = pixel_mm, pixel_mm
        # ncols is nfast, nrows is nslow
        image_size = (
            int(self.header_dict["NCOLS"].split()[0]),
            int(self.header_dict["NROWS"].split()[0]),
        )

        gain = self._calculate_gain(float(self.header_dict["WAVELEN"].split()[0]))
        detector = self._detector_factory.complex(
            "PAD",
            origin.elems,
            fast.elems,
            slow.elems,
            pixel_size,
            image_size,
            (min_trusted_value, full_scale),
        )

        # Here we set specifics, notably parallax correction and
        # QE correction are effectively disabled by setting the simple
        # pixel-to-millimetre strategy and a very high mu value.
        for panel in detector:
            panel.set_gain(gain)
            panel.set_thickness(0.450)  # Assume same as Eiger
            panel.set_material("Si")
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(1e10)

        return detector

    def _beam(self):
        """Make unpolarized beam"""
        wavelength = float(self.header_dict["WAVELEN"].split()[0])

        return self._beam_factory.make_polarized_beam(
            sample_to_source=(0.0, 0.0, 1.0),
            wavelength=wavelength,
            polarization=(0, 1, 0),
            polarization_fraction=0.5,
            probe=Probe.electron,
        )

    def _scan(self):
        start = float(self.header_dict["START"].split()[0])
        incr = float(self.header_dict["INCREME"].split()[0])
        if incr < 0:
            start *= -1
            incr *= -1

        return self._scan_factory.single_file(
            filename=self._image_file,
            exposure_times=1,
            osc_start=start,
            osc_width=incr,
            epoch=None,
        )

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        # It is better to catch FORMAT 86 here and fail with a sensible error msg
        # as soon as something is attempted with the image data rather than in
        # the understand method. Otherwise the user gets FormatBruker reading the
        # image improperly but without failing
        if self.header_dict["FORMAT"] != "100":
            raise NotImplementedError("Only FORMAT 100 images are currently supported")

        f = self.open_file(self._image_file, "rb")
        header_size = int(self.header_dict["HDRBLKS"]) * 512
        f.read(header_size)

        if is_big_endian():
            read_2b = read_uint16_bs
            read_4b = read_uint32_bs
        else:
            read_2b = read_uint16
            read_4b = read_uint32

        # NPIXELB stores the number of bytes/pixel for the data and the underflow
        # table. We expect 1 byte for underflows and either 2 or 1 byte per pixel
        # for the data
        npixelb = [int(e) for e in self.header_dict["NPIXELB"].split()]
        assert npixelb[1] == 1

        if npixelb[0] == 1:
            read_data = read_uint8
        elif npixelb[0] == 2:
            read_data = read_2b
        else:
            raise IncorrectFormatError(f"{npixelb[0]} bytes per pixel is not supported")

        nrows = int(self.header_dict["NROWS"].split()[0])
        ncols = int(self.header_dict["NCOLS"].split()[0])

        raw_data = read_data(streambuf(f), nrows * ncols)

        image_size = (nrows, ncols)
        raw_data.reshape(flex.grid(*image_size))

        (num_underflows, num_2b_overflows, num_4b_overflows) = (
            int(e) for e in self.header_dict["NOVERFL"].split()
        )

        # read underflows
        if num_underflows > 0:
            # stored values are padded to a multiple of 16 bytes
            nbytes = num_underflows + 15 & ~(15)
            underflow_vals = read_uint8(streambuf(f), nbytes)[:num_underflows]
        else:
            underflow_vals = None

        # handle 2 byte overflows
        if num_2b_overflows > 0:
            # stored values are padded to a multiple of 16 bytes
            nbytes = num_2b_overflows * 2 + 15 & ~(15)
            overflow_vals = read_2b(streambuf(f), nbytes // 2)[:num_2b_overflows]
            overflow = flex.int(nrows * ncols, 0)
            sel = (raw_data == 255).as_1d()
            overflow.set_selected(sel, overflow_vals - 255)
            overflow.reshape(flex.grid(*image_size))
            raw_data += overflow

        # handle 4 byte overflows
        if num_4b_overflows > 0:
            # stored values are padded to a multiple of 16 bytes
            nbytes = num_4b_overflows * 4 + 15 & ~(15)
            overflow_vals = read_4b(streambuf(f), nbytes // 4)[:num_4b_overflows]
            overflow = flex.int(nrows * ncols, 0)
            sel = (raw_data == 65535).as_1d()
            overflow.set_selected(sel, overflow_vals - 65535)
            overflow.reshape(flex.grid(*image_size))
            raw_data += overflow

        # handle underflows
        if underflow_vals is not None:
            sel = (raw_data == 0).as_1d()
            underflow = flex.int(nrows * ncols, 0)
            underflow.set_selected(sel, underflow_vals)
            underflow.reshape(flex.grid(*image_size))
            raw_data += underflow

        # handle baseline. num_underflows == -1 means no baseline subtraction. See
        # https://github.com/cctbx/cctbx_project/files/1262952/BISFrameFileFormats.zip
        if num_underflows != -1:
            num_exposures = [int(e) for e in self.header_dict["NEXP"].split()]
            baseline = num_exposures[2]
            raw_data += baseline

        return raw_data
