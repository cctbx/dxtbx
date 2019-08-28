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
