"""
An implementation of the SMV image reader for ADSC images. Inherits from
FormatSMVADSCSN, customised for example on Australian Synchrotron SN 928
which has reversed phi.
"""

import sys

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN


class FormatSMVADSCSN928(FormatSMVADSCSN):
    """A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this, for instrument number 928."""

    @staticmethod
    def understand(image_file):
        """Check to see if this is ADSC SN 928."""

        # check this is detector serial number 928

        size, header = FormatSMVADSCSN.get_smv_header(image_file)
        if int(header["DETECTOR_SN"]) != 928:
            return False

        return True

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header."""

        return self._goniometer_factory.single_axis_reverse()

    def _scan(self):
        """Return the scan information for this image. There may be
        no timestamps in there..."""

        exposure_time = float(self._header_dictionary["TIME"])
        epoch = 0
        osc_start = float(self._header_dictionary["OSC_START"])
        osc_range = float(self._header_dictionary["OSC_RANGE"])

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVADSCSN928.understand(arg))
        print(FormatSMVADSCSN928(arg).get_scan())
        print(FormatSMVADSCSN928(arg).get_goniometer())
        print(FormatSMVADSCSN928(arg).get_detector())
        print(FormatSMVADSCSN928(arg).get_beam())
