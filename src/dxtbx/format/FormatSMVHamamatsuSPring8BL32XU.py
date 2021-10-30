from iotbx.detectors.hamamatsu import HamamatsuImage

from dxtbx.format.FormatSMVHamamatsu import FormatSMVHamamatsu


class FormatSMVHamamatsuSPring8BL32XU(FormatSMVHamamatsu):
    @staticmethod
    def understand(image_file):

        size, header = FormatSMVHamamatsu.get_smv_header(image_file)

        wanted_header_items = ["DETECTOR_NAME"]

        for header_item in wanted_header_items:
            if header_item not in header:
                return 0

        return header["DETECTOR_NAME"] == "Hamamatsu C10158DK"

    def _start(self):

        FormatSMVHamamatsu._start(self)

    def detectorbase_start(self):
        self.detectorbase = HamamatsuImage(self._image_file)
        self.detectorbase.open_file = self.open_file
        self.detectorbase.readHeader()

    def _goniometer(self):
        """Return a model for a simple single-axis reversed direction goniometer."""

        return self._goniometer_factory.single_axis_reverse()
