"""
An implementation of the SMV image reader for ADSC images. Inherits from
FormatSMVADSC, customised for old detector on Diamond Light Source I03,
correctly accounting for the image pedestal & similar
"""

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
