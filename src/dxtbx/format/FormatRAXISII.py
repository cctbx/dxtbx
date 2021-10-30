import sys

from iotbx.detectors.raxis_nonsquare import NonSquareRAXISImage

from dxtbx.format.Format import Format
from dxtbx.format.FormatRAXIS import RAXISHelper


class FormatRAXISII(RAXISHelper, Format):
    @staticmethod
    def understand(image_file):
        try:
            with FormatRAXISII.open_file(image_file, "rb") as fh:
                return fh.read(7) == b"R-AXIS2"
        except OSError:
            return False

    def detectorbase_start(self):
        pass  # override with an empty function

    def _start(self):
        self.detectorbase = NonSquareRAXISImage(self._image_file)
        self.detectorbase.readHeader()

    def _goniometer(self):
        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector"""

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
            exposure_times=1,
            osc_start=self.detectorbase.parameters["OSC_START"],
            osc_width=self.detectorbase.parameters["OSC_RANGE"],
            epoch=None,
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatRAXISII.understand(arg))
