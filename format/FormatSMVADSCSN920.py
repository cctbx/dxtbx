#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for ADSC images. Inherits from
# FormatSMVADSC, customised for old detector on Diamond Light Source I03,
# correctly accounting for the image pedestal & similar

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN


class FormatSMVADSCSN920(FormatSMVADSCSN):
    """A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this, for instrument number 920."""

    @staticmethod
    def understand(image_file):
        """Check to see if this is ADSC SN 920."""

        # check this is detector serial number 920

        size, header = FormatSMVADSCSN.get_smv_header(image_file)

        return int(header["DETECTOR_SN"]) == 920

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatSMVADSCSN.__init__(self, image_file, **kwargs)

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        from boost.python import streambuf
        from dxtbx import read_uint16, read_uint16_bs, is_big_endian
        from scitbx.array_family import flex

        assert len(self.get_detector()) == 1
        image_pedestal = int(self._header_dictionary["IMAGE_PEDESTAL"])
        panel = self.get_detector()[0]
        size = panel.get_image_size()
        if self._header_dictionary["BYTE_ORDER"] == "big_endian":
            big_endian = True
        else:
            big_endian = False

        with FormatSMVADSCSN.open_file(self._image_file, "rb") as fh:
            fh.seek(self._header_size)
            if big_endian == is_big_endian():
                raw_data = read_uint16(streambuf(fh), int(size[0] * size[1]))
            else:
                raw_data = read_uint16_bs(streambuf(fh), int(size[0] * size[1]))

        # apply image pedestal, will result in *negative pixel values*

        raw_data -= image_pedestal

        image_size = panel.get_image_size()
        raw_data.reshape(flex.grid(image_size[1], image_size[0]))

        return raw_data
