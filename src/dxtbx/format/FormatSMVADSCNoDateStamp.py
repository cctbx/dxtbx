import os
import sys

from iotbx.detectors.adsc import ADSCImage

from dxtbx.format.FormatSMVADSC import FormatSMVADSC


class FormatSMVADSCNoDateStamp(FormatSMVADSC):
    """A class for reading SMV format ADSC images, with detector serial number"""

    @staticmethod
    def understand(image_file):

        # assert for this that the image file has to be a file not a URL
        if not os.path.exists(image_file):
            return False

        size, header = FormatSMVADSC.get_smv_header(image_file)

        wanted_header_items = ["TIME"]
        if any(item not in header for item in wanted_header_items):
            return False

        unwanted_header_items = ["DATE"]

        if any(item in header for item in unwanted_header_items):
            return False

        return True

    def _scan(self):
        """Return the scan information for this image, using the timestamp
        from the file rather than the header."""

        exposure_time = float(self._header_dictionary["TIME"])

        epoch = float(os.stat(self._image_file)[8])

        osc_start = float(self._header_dictionary["OSC_START"])
        osc_range = float(self._header_dictionary["OSC_RANGE"])

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )

    def detectorbase_start(self):
        self.detectorbase = ADSCImage(self._image_file)
        self.detectorbase.open_file = self.open_file
        self.detectorbase.readHeader()


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVADSCNoDateStamp.understand(arg))
