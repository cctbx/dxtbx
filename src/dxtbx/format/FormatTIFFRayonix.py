"""An implementation of the TIFF image reader for Rayonix images."""


import datetime
import struct
import sys
import time

from boost_adaptbx.boost.python import streambuf
from iotbx.detectors.mar import MARImage
from scitbx.array_family import flex

from dxtbx.ext import read_uint16
from dxtbx.format.FormatTIFF import FormatTIFF


class FormatTIFFRayonix(FormatTIFF):
    """A class for reading TIFF format Rayonix images, and correctly
    constructing a model for the experiment from this."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Rayonix TIFF format image,
        i.e. we can make sense of it. This simply checks that records which
        describe the size of the image match with the TIFF records which do
        the same."""

        width, height, depth, order, bytes = FormatTIFF.get_tiff_header(image_file)

        if len(bytes) != 4096:
            return False

        if order == FormatTIFF.LITTLE_ENDIAN:
            endian = "<"
        else:
            endian = ">"

        _I = endian + "I"

        _width = struct.unpack(_I, bytes[1024 + 80 : 1024 + 84])[0]
        _height = struct.unpack(_I, bytes[1024 + 84 : 1024 + 88])[0]
        _depth = struct.unpack(_I, bytes[1024 + 88 : 1024 + 92])[0]

        if width != _width or height != _height or depth != _depth:
            return False

        # pretty sure all MARCCD / Rayonix detectors are square...
        if width != height:
            return False

        nimages = struct.unpack(_I, bytes[1024 + 112 : 1024 + 116])[0]
        origin = struct.unpack(_I, bytes[1024 + 116 : 1024 + 120])[0]
        orientation = struct.unpack(_I, bytes[1024 + 120 : 1024 + 124])[0]
        view = struct.unpack(_I, bytes[1024 + 124 : 1024 + 128])[0]

        if nimages != 1 or origin != 0 or orientation != 0 or view != 0:
            return False

        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        width, height, depth, order, bytes = FormatTIFF.get_tiff_header(image_file)

        # comment block - where the detector serial number may (or may not) be stored
        # comments = bytes[1024+1440:1024+1440+512]

        self._header_size = 4096

        if order == FormatTIFF.LITTLE_ENDIAN:
            self._I = "<I"
            self._i = "<i"
            self._ii = "<ii"
        else:
            self._I = ">I"
            self._i = ">i"
            self._ii = ">ii"

        super().__init__(image_file, **kwargs)

    def detectorbase_start(self):
        self.detectorbase = MARImage(self._image_file)
        self.detectorbase.readHeader()

    # FIXME have implemented none of those which follow...

    def _goniometer(self):
        """Return a model for goniometer corresponding to the values stored
        in the image header. In the first instance assume this is a single
        axis and raise exception otherwise."""

        starts, ends, offset, width = self._get_rayonix_scan_angles()

        for j, starts_j in enumerate(starts):
            if j != offset:
                assert starts_j == ends[j]

        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector, which at the moment [NKS does not] insist
        that the offsets and rotations are all 0.0."""

        starts, ends, offset, width = self._get_rayonix_scan_angles()
        rotations = self._get_rayonix_detector_rotations()

        # NKS removed the assertion that two-theta offset is 0.0; support the general case
        assert starts[0] == ends[0]
        two_theta = starts[0]

        # assert that the rotations are all 0.0

        assert rotations[0] == 0.0
        assert rotations[1] == 0.0
        assert rotations[2] == 0.0

        distance = self._get_rayonix_distance()
        beam_x, beam_y = self._get_rayonix_beam_xy()
        pixel_size = self._get_rayonix_pixel_size()
        image_size = self._tiff_width, self._tiff_height
        overload = struct.unpack(self._i, self._tiff_header_bytes[1128:1132])[0]
        underload = 0

        bias = int(round(self._get_rayonix_bias()))
        underload -= bias
        overload -= bias

        beam = beam_x * pixel_size[0], beam_y * pixel_size[1]

        return self._detector_factory.two_theta(
            "CCD",
            distance,
            beam,
            "+x",
            "-y",
            "+x",
            two_theta,
            pixel_size,
            image_size,
            (underload, overload),
            [],
            gain=self._get_rayonix_gain(),
        )

    def _beam(self):
        """Return a simple model for the beam."""

        wavelength = (
            struct.unpack(self._i, self._tiff_header_bytes[1932:1936])[0] * 1.0e-5
        )

        return self._beam_factory.simple(wavelength)

    def _scan(self):
        """Return the scan information for this image."""

        exposure_time = self._get_rayonix_times()[1]
        epoch = time.mktime(self._get_rayonix_timestamp())

        starts, ends, offset, width = self._get_rayonix_scan_angles()

        osc_start = starts[offset]
        osc_range = width

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )

    ####################################################################
    #                                                                  #
    # Helper methods to get all of the values out of the TIFF header   #
    # - separated out to assist with code clarity                      #
    #                                                                  #
    ####################################################################

    def _get_rayonix_distance(self):
        """Look in the usual places for the detector distance, return this
        as a float in mm."""

        distance = struct.unpack(self._i, self._tiff_header_bytes[1664:1668])[0]

        if distance != 0:
            return distance * 0.001

        distance = struct.unpack(self._i, self._tiff_header_bytes[1720:1724])[0]

        if distance != 0:
            return distance * 0.001

        raise RuntimeError("cannot find distance in header")

    def _get_rayonix_beam_xy(self):
        """Get the beam x, y positions which are defined in the standard
        to be in pixels. X and Y are not defined by the documentation, so
        taking as a given that these are horizontal and vertical. N.B.
        the documentation states that the horizontal direction is fast."""

        beam_x, beam_y = struct.unpack(self._ii, self._tiff_header_bytes[1668:1676])[:2]

        return beam_x * 0.001, beam_y * 0.001

    def _get_rayonix_pixel_size(self):
        """Get the pixel sizes in mm."""

        pixel_x, pixel_y = struct.unpack(self._ii, self._tiff_header_bytes[1796:1804])[
            :2
        ]

        return pixel_x * 1.0e-6, pixel_y * 1.0e-6

    def _get_rayonix_times(self):
        """Get the integration, exposure times in seconds."""

        integration, exposure = struct.unpack(
            self._ii, self._tiff_header_bytes[1676:1684]
        )[:2]

        return integration * 0.001, exposure * 0.001

    def _get_rayonix_timestamp(self):
        """Get the image acquisition timestamp."""

        timestamp = self._tiff_header_bytes[2368:2400]

        month = int(timestamp[:2])
        day = int(timestamp[2:4])
        hour = int(timestamp[4:6])
        minute = int(timestamp[6:8])
        year = int(timestamp[8:12])
        second = int(timestamp[13:15])
        return datetime.datetime(year, month, day, hour, minute, second).timetuple()

    def _get_rayonix_scan_angles(self):
        """Get the scan angles for: twotheta, omega, chi, kappa, phi, delta,
        gamma. The exact definitions for these are somewhat poorly defined,
        though I presume that they will come from some kind of standard
        goniometer definition... Also returns the scan axis offset and
        the apparent oscillation width, as these are sometimes not properly
        recorded elsewhere."""

        if self._tiff_byte_order == FormatTIFF.LITTLE_ENDIAN:
            iiiiiii = "<iiiiiii"
        else:
            iiiiiii = ">iiiiiii"

        start_angles = struct.unpack(iiiiiii, self._tiff_header_bytes[1692:1720])[:7]

        end_angles = struct.unpack(iiiiiii, self._tiff_header_bytes[1724:1752])[:7]

        axis_offset = struct.unpack(self._i, self._tiff_header_bytes[1756:1760])[0]

        axis_range = struct.unpack(self._i, self._tiff_header_bytes[1760:1764])[0]

        starts_degrees = [s * 0.001 for s in start_angles]
        ends_degrees = [e * 0.001 for e in end_angles]

        return starts_degrees, ends_degrees, axis_offset, axis_range * 0.001

    def _get_rayonix_detector_rotations(self):
        """Get the recorded rotx, roty, rotz of the detector - which in most
        cases will probably all be 0.0."""

        if self._tiff_byte_order == FormatTIFF.LITTLE_ENDIAN:
            iii = "<iii"
        else:
            iii = ">iii"

        rot_angles = struct.unpack(iii, self._tiff_header_bytes[1764:1776])[:3]

        rot_degrees = [r * 0.001 for r in rot_angles]

        return rot_degrees

    def _get_rayonix_bias(self):
        """Get the image bias."""

        bias = struct.unpack(self._i, self._tiff_header_bytes[1804:1808])[0]

        return bias * 1.0e-3

    def _get_rayonix_gain(self, model=None):
        """Return an appropriate gain value in ADU per captured X-ray for a
        Rayonix CCD module. If the model is None (unknown) then make a guess based
        on header details"""

        # Values are based on a combination of manufacturer datasheets,
        # James Holton's list at http://bl831.als.lbl.gov/xtalsize.html and
        # empirical guesswork based on spot finding results. The accuracy should
        # not be expected to be better than 20% and indeed may be worse than that.
        # Values may refer to gain in ADU per incident photon rather than the
        # more appropriate ADU per captured photon.

        pixel_size = self._get_rayonix_pixel_size()
        image_size = self._tiff_width, self._tiff_height

        if model is None:
            model = "UNKNOWN"
            # Guesswork based on information from the header
            if image_size == (2048, 2048) and pixel_size == (0.079348, 0.079348):
                model = "165"  # 2x2 binning
            elif image_size == (3072, 3072):
                model = "225"  # 2x2 binning, 73 micron pixels. But might be 225HE?
            elif image_size == (4096, 4096):
                if pixel_size == (0.079346, 0.079346):
                    model = "325"  # 2x2 binning, 79 micron pixels. But might be 325HE?
                elif pixel_size == (0.073242, 0.073242):
                    model = "300"  # 2x2 binning, 73 micron pixels. But might be 300HE?

        model = model.upper()
        if model.startswith("165"):
            return 1.0
        if model.startswith("225HE"):
            return 2.8  # This is inferred based on the value for 300HE
        if model.startswith("225"):
            return 1.5
        if model.startswith("325HE"):
            return 2.4
        if model.startswith("325"):
            return 1.3
        if model.startswith("300HE"):
            return 2.8
        if model.startswith("300"):
            return 1.5

        # Unidentified model
        return 1.0

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        # currently have no non-little-endian machines...

        assert self._tiff_byte_order == FormatTIFF.LITTLE_ENDIAN

        bias = int(round(self._get_rayonix_bias()))

        assert len(self.get_detector()) == 1
        size = self.get_detector()[0].get_image_size()
        f = FormatTIFF.open_file(self._image_file)
        f.read(self._header_size)
        raw_data = read_uint16(streambuf(f), int(size[0] * size[1])) - bias
        image_size = self.get_detector()[0].get_image_size()
        raw_data.reshape(flex.grid(image_size[1], image_size[0]))

        return raw_data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatTIFFRayonix.understand(arg))
