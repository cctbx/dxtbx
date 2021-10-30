"""
Implementation of an ImageFormat class to read TIFF format image but not -
in the first instance - actually provide a full image representation. This
is simply there to set everything up for the Mar / Rayonix CCD readers
which really will acquire the full image including header information
and generate the experimental model representations.
"""


from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format
from dxtbx.format.FormatTIFFHelpers import (
    BIG_ENDIAN,
    LITTLE_ENDIAN,
    read_basic_tiff_header,
)


class FormatTIFF(Format):
    """An image reading class for TIFF format images i.e. those from Dectris
    and Rayonix, which start with a standard TIFF header (which is what is
    handled here) and have their own custom header following, which must
    be handled by the inheriting subclasses."""

    LITTLE_ENDIAN = LITTLE_ENDIAN
    BIG_ENDIAN = BIG_ENDIAN

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an TIFF format image, i.e. we can
        make sense of it."""

        try:
            return bool(read_basic_tiff_header(image_file))
        except Exception:
            return False

    @staticmethod
    def get_tiff_header(image_file):
        """Pun to get to the image header etc."""

        width, height, depth, header, order = read_basic_tiff_header(image_file)

        header_bytes = FormatTIFF.open_file(image_file, "rb").read(header)

        return width, height, depth // 8, order, header_bytes

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        super().__init__(image_file, **kwargs)

    def _start(self):
        """Open the image file, read the image header, copy it into memory
        for future inspection."""

        width, height, depth, header, order = read_basic_tiff_header(self._image_file)

        self._tiff_width = width
        self._tiff_height = height
        self._tiff_depth = depth // 8
        self._tiff_header_bytes = FormatTIFF.open_file(self._image_file, "rb").read(
            header
        )
        self._tiff_byte_order = order
