"""
Sub class of FormatTIFFRayonix specialized for the XPP Rayonix dectector at LCLS

Images from the XPP Rayonix detector have several unitialized values, such as
distance, wavelength, etc.  Set these values to zero so the images can be at
least viewed.
"""


import re
import struct
import sys

from iotbx.detectors.mar import MARImage

from dxtbx.format.FormatTIFFRayonix import FormatTIFFRayonix


def check(l):
    """Sets l or values in l that are less than zero to zero"""
    if not isinstance(l, list) and not isinstance(l, tuple):
        if l < 0:
            return 0
        else:
            return l
    ret = []
    for val in l:
        if isinstance(val, list):
            ret.append(check(val))
        else:
            if val < 0:
                val = 0
            ret.append(val)

    if isinstance(l, tuple):
        return tuple(ret)
    return ret


_xpp_pattern = re.compile(r".*xpp[a-zA-Z][0-9][0-9][0-9][0-9].*")


class FormatTIFFRayonixXPP(FormatTIFFRayonix):
    """A class for reading TIFF format Rayonix images, and correctly
    constructing a model for the experiment from this."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an XPP Rayonix TIFF"""

        with open(image_file, "rb") as fh:

            def get(index, typelen, offset=1024):
                fh.seek(offset + index)
                rawdata = fh.read(typelen)
                return struct.unpack("<c", rawdata)[0]

            data = [get(1152 + i, 1) for i in range(128)]
        filepath = b"".join(c for c in data if 32 < ord(c) < 127).decode(
            "latin-1", "replace"
        )

        return bool(_xpp_pattern.match(filepath))

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        MARImage._read_header_asserts = lambda self: None
        super().__init__(image_file, **kwargs)

    ####################################################################
    #                                                                  #
    # Helper methods to get all of the values out of the TIFF header   #
    # - separated out to assist with code clarity                      #
    #                                                                  #
    ####################################################################

    def _get_rayonix_beam_xy(self):
        """Get the beam x, y positions which are defined in the standard
        to be in pixels. X and Y are not defined by the documentation, so
        taking as a given that these are horizontal and vertical. N.B.
        the documentation states that the horizontal direction is fast."""

        beam_x, beam_y = struct.unpack(self._ii, self._tiff_header_bytes[1668:1676])[:2]
        pixel_x, pixel_y = struct.unpack(self._ii, self._tiff_header_bytes[1796:1804])[
            :2
        ]

        return beam_x * 1000 / pixel_x, beam_y * 1000 / pixel_y

    def _get_rayonix_detector_rotations(self):
        return check(FormatTIFFRayonix._get_rayonix_detector_rotations(self))

    def _get_rayonix_distance(self):
        return check(FormatTIFFRayonix._get_rayonix_distance(self))

    def _get_rayonix_pixel_size(self):
        return check(FormatTIFFRayonix._get_rayonix_pixel_size(self))

    def _get_rayonix_times(self):
        return check(FormatTIFFRayonix._get_rayonix_times(self))

    def _get_rayonix_timestamp(self):
        return check(FormatTIFFRayonix._get_rayonix_timestamp(self))

    def _get_rayonix_scan_angles(self):
        return check(FormatTIFFRayonix._get_rayonix_scan_angles(self))


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatTIFFRayonixXPP.understand(arg))
