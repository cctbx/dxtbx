from iotbx.detectors.hamamatsu import HamamatsuImage

from dxtbx.format.FormatSMVADSC import FormatSMVADSC


class FormatSMVHamamatsu(FormatSMVADSC):
    @staticmethod
    def understand(image_file):

        size, header = FormatSMVHamamatsu.get_smv_header(image_file)

        wanted_header_items = ["DETECTOR_NAME"]
        if any(item not in header for item in wanted_header_items):
            return False

        return "hamamatsu" in header["DETECTOR_NAME"].lower()

    def detectorbase_start(self):
        self.detectorbase = HamamatsuImage(self._image_file)
        self.detectorbase.open_file = self.open_file
        self.detectorbase.readHeader()
