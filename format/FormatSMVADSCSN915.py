#!/usr/bin/env python
# FormatSMVADSCSN915.py
#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# ADSC SMV Format for Q315 SN 915, installed at BL38B1 at SPring-8. Resembles
# but FormatSMVADSCSN920 but returns a reverse phi goniometer

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN


class FormatSMVADSCSN915(FormatSMVADSCSN):
    """A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this, for instrument number 915."""

    @staticmethod
    def understand(image_file):
        """Check to see if this is ADSC SN 915."""

        # check this is detector serial number 915

        size, header = FormatSMVADSCSN.get_smv_header(image_file)

        if int(header["DETECTOR_SN"]) != 915:
            return False

        return True

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
        f = FormatSMVADSCSN.open_file(self._image_file, "rb")
        f.read(self._header_size)

        if self._header_dictionary["BYTE_ORDER"] == "big_endian":
            big_endian = True
        else:
            big_endian = False

        if big_endian == is_big_endian():
            raw_data = read_uint16(streambuf(f), int(size[0] * size[1]))
        else:
            raw_data = read_uint16_bs(streambuf(f), int(size[0] * size[1]))

        # apply image pedestal, will result in *negative pixel values*

        raw_data -= image_pedestal

        image_size = panel.get_image_size()
        raw_data.reshape(flex.grid(image_size[1], image_size[0]))

        return raw_data

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer with reversed
        direction."""

        return self._goniometer_factory.single_axis_reverse()


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatSMVADSC.understand(arg))
