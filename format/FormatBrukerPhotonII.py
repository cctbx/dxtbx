import sys

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


class FormatBrukerPhotonII(FormatBruker):
    @staticmethod
    def understand(image_file):

        try:
            header_lines = FormatBruker.read_header_lines(image_file)
        except OSError:
            return False

        header_dic = FormatBruker.parse_header(header_lines)

        dettype = header_dic.get("DETTYPE")
        if dettype is None:
            return False
        if not dettype.startswith("CMOS-PHOTONII"):
            return False

        return True

    def _start(self):

        try:
            header_lines = FormatBruker.read_header_lines(self._image_file)
        except OSError:
            return False

        self.header_dict = FormatBrukerPhotonII.parse_header(header_lines)

        # The Photon II format can't currently use BrukerImage, see
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
            scan_axis = "OMEGA"  # default

        # https://journals.iucr.org/d/issues/2014/10/00/dz5309/dz5309sup1.pdf
        axes = flex.vec3_double(((0, -1, 0), (0, 0, 1), (0, 1, 0)))
        omega -= 180
        angles = flex.double((phi, chi, omega))

        return self._goniometer_factory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis
        )

    @staticmethod
    def _estimate_gain(wavelength):
        """Estimate the detector gain based on values provided by Bruker. Each ADU
        corresponds to 36.6 electrons. The X-ray conversion results in deposited
        charge according to the following table for typical home sources:

        In (0.5136 A): 359.6893 e/X-ray
        Ag (0.5609 A): 329.3748 e/X-ray
        Mo (0.7107 A): 259.9139 e/X-ray
        Ga (1.3414 A): 137.6781 e/X-ray
        Cu (1.5418 A): 119.8156 e/X-ray

        This fits the linear model (1/G) = -0.0000193358 + 0.1981607255 * wavelength
        extremely well.
        """
        inv_gain = -0.0000193358 + 0.1981607255 * wavelength
        assert inv_gain > 0.1
        return 1.0 / inv_gain

    def _detector(self):
        # goniometer angles in ANGLES are 2-theta, omega, phi, chi (FIXED)
        two_theta = float(self.header_dict["ANGLES"].split()[0])

        overload = float(self.header_dict["CCDPARM"].split()[-1])
        underload = -1

        fast = matrix.col((1, 0, 0))
        slow = matrix.col((0, 1, 0))
        beam = matrix.col((0, 0, 1))
        pixel_mm = 5.0 / float(self.header_dict["DETTYPE"].split()[1])
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

        # Not a CCD, but is an integrating detector. Photon II has a 90 um Gadox
        # scintillator.
        gain = self._estimate_gain(float(self.header_dict["WAVELEN"].split()[0]))
        return self._detector_factory.complex(
            "CCD",
            origin.elems,
            fast.elems,
            slow.elems,
            pixel_size,
            image_size,
            (underload, overload),
            gain=gain,
        )

    def _beam(self):
        wavelength = float(self.header_dict["WAVELEN"].split()[0])

        return self._beam_factory.simple(wavelength)

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
            raise NotImplementedError(
                "Only FORMAT 100 images from the Photon II are currently supported"
            )

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
            raise IncorrectFormatError(
                "{} bytes per pixel is not supported".format(npixelb[0])
            )

        nrows = int(self.header_dict["NROWS"].split()[0])
        ncols = int(self.header_dict["NCOLS"].split()[0])

        raw_data = read_data(streambuf(f), nrows * ncols)

        image_size = (nrows, ncols)
        raw_data.reshape(flex.grid(*image_size))

        (num_underflows, num_2b_overflows, num_4b_overflows) = [
            int(e) for e in self.header_dict["NOVERFL"].split()
        ]

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


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatBrukerPhotonII.understand(arg))
