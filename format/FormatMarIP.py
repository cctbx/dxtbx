"""An ImageFormat class to read MarIP-format image"""

from __future__ import absolute_import, division, print_function

import sys

from iotbx.detectors.detectorbase import DetectorImageBase

from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format

try:
    # Try to import the modernised pycbf
    import pycbf
    from pycbf.img import Img

    import scitbx.array_family.flex as flex

    print(f"Using modernised pycbf-{pycbf.__version__}")
except ModuleNotFoundError:
    Img = None
    # Replicate fallback logic from iotbx.detectors.marIP
    try:
        from cbflib_adaptbx import Mar345Adaptor
    except ModuleNotFoundError:
        from iotbx.detectors.marIP import NullAdaptor as Mar345Adaptor

if Img is not None:
    # We have the modernised pycbf with the img bindings.
    # Declare Mar345Adaptor in pure python
    class Mar345Adaptor:  # noqa: F811
        def __init__(self, filename):
            self.filename = filename

            self._read_header = False
            self._read_data = False
            self._data_position = None
            self._metadata = None
            self._img = Img()
            self._img.set_tags(0)
            self._img.set_dimensions(0, 0)

        def read_header(self):
            if self._read_header:
                return
            with open(self.filename, "rb") as f:
                self._metadata = self._img.read_mar345header(f)
                self._data_position = f.tell()
            detector = self._img.get_field("DETECTOR")
            if "mar" not in detector.lower() or "345" not in detector.lower():
                raise RuntimeError(
                    f"Detector type other than mar345 from mar345 reader ({detector})"
                )

        def read_data(self):
            if self._read_data:
                return
            self.read_header()
            with open(self.filename, "rb") as f:
                f.seek(self._data_position)
                self._img.read_mar345data(f, self._metadata)

        def rawdata(self):
            self.read_data()
            return flex.int(self._img.image)

        def columns(self):
            return self._img.columns

        def rows(self):
            return self._img.rows

        def size1(self):
            self.read_header()
            return self._metadata[0]

        def size2(self):
            self.read_header()
            return self._metadata[2]

        def pixel_size(self):
            self.read_header()
            return self._img.get_number("PIXEL SIZE")

        def wavelength(self):
            self.read_header()
            return self._img.get_number("WAVELENGTH")

        def distance(self):
            self.read_header()
            return self._img.get_number("DISTANCE")

        def gain(self):
            return 1.55

        def overload(self):
            # Determined from a single image.
            return 249862

        def osc_range(self):
            self.read_header()
            return self._img.get_number("OSCILLATION RANGE")

        def osc_start(self):
            self.read_header()
            return self._img.get_number("PHI")

        def twotheta(self):
            self.read_header()
            return self._img.get_number("TWOTHETA")

        def exposure_time(self):
            self.read_header()
            return self._img.get_number("EXPOSURE_TIME")


class MARIPImage(DetectorImageBase):
    def __init__(self, filename):
        DetectorImageBase.__init__(self, filename)
        self.adaptor = Mar345Adaptor(filename)
        self.vendortype = "MARIP"

    def readHeader(self):
        self.parameters = {
            "SIZE1": self.adaptor.size1(),
            "SIZE2": self.adaptor.size2(),
            "CCD_IMAGE_SATURATION": self.adaptor.overload(),
            "PIXEL_SIZE": self.adaptor.pixel_size(),
            "OSC_START": self.adaptor.osc_start(),
            "DISTANCE": self.adaptor.distance(),
            "WAVELENGTH": self.adaptor.wavelength(),
            "BEAM_CENTER_X": self.beam_center_slow(),
            "BEAM_CENTER_Y": self.beam_center_fast(),
            "OSC_RANGE": self.adaptor.osc_range(),
            "TWOTHETA": self.adaptor.twotheta(),
            "DETECTOR_SN": 0,
        }

    def beam_center_slow(self):
        return self.adaptor.size1() * self.adaptor.pixel_size() / 2.0

    def beam_center_fast(self):
        return self.adaptor.size2() * self.adaptor.pixel_size() / 2.0

    def fileLength(self):
        return 0

    def getEndian(self):
        return 0

    def read(self):
        self.bin_safe_set_data(self.adaptor.rawdata())

    def dataoffset(self):
        return 0

    def integerdepth(self):
        return 0


class FormatMarIP(Format):
    """An image reading class for MarIP-format images
    Positive identification:  first 140 bytes contain the string "mar research"
    """

    @staticmethod
    def understand(image_file):
        try:
            with FormatMarIP.open_file(image_file, "rb") as fh:
                return b"mar research" in fh.read(140)
        except IOError:
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
