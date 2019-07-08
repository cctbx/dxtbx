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

        return int(header["DETECTOR_SN"]) == 915

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

        assert len(self.get_detector()) == 1
        panel = self.get_detector()[0]
        image_size = panel.get_image_size()
        raw_data = self._get_endianic_raw_data(size=image_size)

        # apply image pedestal, will result in *negative pixel values*
        image_pedestal = int(self._header_dictionary["IMAGE_PEDESTAL"])
        raw_data -= image_pedestal

        return raw_data

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer with reversed
        direction."""

        return self._goniometer_factory.single_axis_reverse()
