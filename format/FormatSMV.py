"""
Implementation of an ImageFormat class to read SMV format image

but not - in the first instance - actually provide a full image
representation. This is simply there to set everything up for the ADSC
and Rigaku Saturn image readers which really will acquire the full image
including header information and generate the experimental model
representations.
"""


from boost_adaptbx.boost.python import streambuf
from scitbx.array_family import flex

from dxtbx import IncorrectFormatError
from dxtbx.ext import is_big_endian, read_uint16, read_uint16_bs
from dxtbx.format.Format import Format


class FormatSMV(Format):
    """An image reading class for SMV format images i.e. those from ADSC and
    Rigaku which start with:

    {
    HEADER_BYTES=  512;

    and contain a list of keyword-value pairs thereafter which define the
    header. The keywords etc. for these will depend on the instrument
    manufacturer and will be interpreted by subclasses of this class. Note
    also that every line is finished with a semicolon."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an SMV format image, i.e. we can
        make sense of it."""
        with FormatSMV.open_file(image_file, "rb") as fh:
            return fh.read(15) == b"{\nHEADER_BYTES="

    @staticmethod
    def get_smv_header(image_file):
        with FormatSMV.open_file(image_file, "rb") as fh:
            header_info = fh.read(45).decode("ascii", "ignore")
            header_size = int(
                header_info.split("\n")[1].split("=")[1].replace(";", "").strip()
            )
            fh.seek(0)
            header_text = fh.read(header_size).decode("ascii", "ignore")
        header_dictionary = {}

        # Check that we have the whole header, contained within { }.  Stop
        # extracting data once a record solely composed of a closing curly
        # brace is seen.  If there is no such character in header_text
        # either HEADER_BYTES caused a short read of the header or the
        # header is malformed.
        for record in header_text.split("\n"):
            if record == "}":
                break
            if "=" not in record:
                continue

            key, value = record.replace(";", "").split("=")

            header_dictionary[key.strip()] = value.strip()

        return header_size, header_dictionary

    def _get_endianic_raw_data(self, size):
        big_endian = self._header_dictionary["BYTE_ORDER"] == "big_endian"

        with self.open_file(self._image_file, "rb") as fh:
            fh.seek(self._header_size)

            if big_endian == is_big_endian():
                raw_data = read_uint16(streambuf(fh), int(size[0] * size[1]))
            else:
                raw_data = read_uint16_bs(streambuf(fh), int(size[0] * size[1]))

        # note that x and y are reversed here
        raw_data.reshape(flex.grid(size[1], size[0]))
        return raw_data

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        super().__init__(image_file, **kwargs)

    def _start(self):
        """Open the image file, read the image header, copy the key / value
        pairs into an internal dictionary self._header_dictionary along with
        the length of the header in bytes self._header_size."""

        self._header_size, self._header_dictionary = FormatSMV.get_smv_header(
            self._image_file
        )

    def get_beam_direction(self):
        return [
            float(sv) for sv in self._header_dictionary["SOURCE_VECTORS"].split()[:3]
        ]

    def get_beam_pixels(self, detector_name):
        return [
            float(bp)
            for bp in self._header_dictionary[
                "%sSPATIAL_DISTORTION_INFO" % detector_name
            ].split()[:2]
        ]

    def get_beam_polarization(self):
        polarization = [
            float(sp) for sp in self._header_dictionary["SOURCE_POLARZ"].split()
        ]
        p_fraction = polarization[0]
        p_plane = polarization[1:]
        return p_fraction, p_plane

    def get_detector_axes(self, detector_name):
        return [
            float(v)
            for v in self._header_dictionary[
                "%sDETECTOR_VECTORS" % detector_name
            ].split()
        ]

    def get_distortion(self, detector_name):
        return [
            float(sdv)
            for sdv in self._header_dictionary[
                "%sSPATIAL_DISTORTION_VECTORS" % detector_name
            ].split()
        ]

    def get_gonio_axes(self, detector_name):
        return [
            float(gv)
            for gv in self._header_dictionary["%sGONIO_VECTORS" % detector_name].split()
        ]

    def get_gonio_values(self, detector_name):
        return [
            float(gv)
            for gv in self._header_dictionary["%sGONIO_VALUES" % detector_name].split()
        ]

    def get_image_size(self, detector_name):
        return [
            int(dd)
            for dd in self._header_dictionary[
                "%sDETECTOR_DIMENSIONS" % detector_name
            ].split()
        ]

    def get_pixel_size(self, detector_name):
        return [
            float(ps)
            for ps in self._header_dictionary[
                "%sSPATIAL_DISTORTION_INFO" % detector_name
            ].split()[2:]
        ]

    def get_rotation(self):
        return [float(r) for r in self._header_dictionary["ROTATION"].split()]
