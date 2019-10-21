"""
ADSC SMV Format for Q315 SN 915, installed at BL38B1 at SPring-8. Resembles
but FormatSMVADSCSN920 but returns a reverse phi goniometer
"""
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

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        assert len(self.get_detector()) == 1
        panel = self.get_detector()[0]
        image_size = panel.get_image_size()
        return self._get_endianic_raw_data(size=image_size)

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer with reversed
        direction."""

        return self._goniometer_factory.single_axis_reverse()
