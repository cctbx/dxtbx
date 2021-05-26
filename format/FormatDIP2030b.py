import sys

from iotbx.detectors.macscience import DIPImage

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format


class FormatDIP2030b(Format):
    @staticmethod
    def understand(image_file):
        # for MacScience DIP2030b only, file size is exactly 18001024 bytes
        headerstart = 3000 * 3000 * 2
        try:
            with FormatDIP2030b.open_file(image_file, "rb") as fh:
                fh.seek(headerstart)
                rawheader = fh.read(1024)
                eof = fh.read(1)  # end of file
        except OSError:
            return False

        return eof == b"" and rawheader[0:3] == b"DIP"

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        super().__init__(image_file, **kwargs)

    def detectorbase_start(self):
        pass

    def _start(self):
        self.detectorbase = DIPImage(self._image_file)
        self.detectorbase.readHeader()

    def _goniometer(self):
        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector"""

        # twotheta = self.detectorbase.parameters["TWOTHETA"]
        # At present, ignore non-zero two theta for the dxtbx model
        # XXX Return to this issue later.
        return self._detector_factory.simple(
            sensor="IMAGE_PLATE",
            distance=self.detectorbase.distance,
            beam_centre=(self.detectorbase.beamx, self.detectorbase.beamy),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(self.detectorbase.pixel_size, self.detectorbase.pixel_size),
            image_size=(self.detectorbase.size1, self.detectorbase.size2),
            trusted_range=(0, self.detectorbase.saturation),
            mask=[],
        )  # a list of dead rectangles

    def _beam(self):
        """Return a simple model for the beam."""

        return self._beam_factory.simple(self.detectorbase.wavelength)

    def _scan(self):
        """Return the scan information for this image."""

        return self._scan_factory.single_file(
            filename=self._image_file,
            exposure_times=self.detectorbase.parameters["TIME"],
            osc_start=self.detectorbase.osc_start,
            osc_width=self.detectorbase.deltaphi,
            epoch=None,
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatDIP2030b.understand(arg))
