"""An ImageFormat class to read MarIP-format image"""


import sys

from iotbx.detectors.marIP import MARIPImage

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format


class FormatMarIP(Format):
    """An image reading class for MarIP-format images
    Positive identification:  first 140 bytes contain the string "mar research"
    """

    @staticmethod
    def understand(image_file):
        try:
            with FormatMarIP.open_file(image_file, "rb") as fh:
                return b"mar research" in fh.read(140)
        except OSError:
            return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        Format.__init__(self, image_file, **kwargs)

    def detectorbase_start(self):
        pass

    def _start(self):
        """Open the image file, read the image header, copy the key / value
        pairs into an internal dictionary self._header_dictionary along with
        the length of the header in bytes self._header_size."""
        self.detectorbase = MARIPImage(self._image_file)
        self.detectorbase.readHeader()

    def _goniometer(self):
        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector, which at the moment insists
        that the offsets and rotations are all 0.0."""

        assert self.detectorbase.parameters["TWOTHETA"] == 0.0

        return self._detector_factory.simple(
            sensor="IMAGE_PLATE",
            distance=self.detectorbase.parameters["DISTANCE"],
            beam_centre=(
                self.detectorbase.parameters["BEAM_CENTER_X"],
                self.detectorbase.parameters["BEAM_CENTER_Y"],
            ),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(
                self.detectorbase.parameters["PIXEL_SIZE"],
                self.detectorbase.parameters["PIXEL_SIZE"],
            ),
            image_size=(
                self.detectorbase.parameters["SIZE1"],
                self.detectorbase.parameters["SIZE2"],
            ),
            trusted_range=(0, self.detectorbase.parameters["CCD_IMAGE_SATURATION"]),
            mask=[],
        )

    def _beam(self):
        """Return a simple model for the beam."""

        return self._beam_factory.simple(self.detectorbase.parameters["WAVELENGTH"])

    def _scan(self):
        """Return the scan information for this image."""

        return self._scan_factory.single_file(
            filename=self._image_file,
            exposure_times=self.detectorbase.adaptor.exposure_time(),
            osc_start=self.detectorbase.parameters["OSC_START"],
            osc_width=self.detectorbase.parameters["OSC_RANGE"],
            epoch=None,
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatMarIP.understand(arg))
