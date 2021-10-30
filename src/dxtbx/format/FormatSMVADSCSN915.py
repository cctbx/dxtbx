"""
ADSC SMV Format for Q315 SN 915, installed at BL38B1 at SPring-8. Resembles
FormatSMVADSCSN920 but returns a reverse phi goniometer
"""

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

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer with reversed
        direction."""

        return self._goniometer_factory.single_axis_reverse()
