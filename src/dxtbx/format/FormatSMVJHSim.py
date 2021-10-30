"""An implementation of the SMV image reader for JHSim images."""


import calendar
import sys
import time

from iotbx.detectors import SMVImage

from dxtbx.format.FormatSMV import FormatSMV


class FormatSMVJHSim(FormatSMV):
    """A class for reading SMV format JHSim images, and correctly constructing
    a model for the experiment from this."""

    # all ADSC detectors generate images with an ADC offset of 40
    # for Mar/Rayonix it is 10
    # Rigaku SMV uses 20, and 5 for image plate formats
    # for one particular simulation, I used 1
    ADC_OFFSET = 1
    image_pedestal = 1

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an JHSim SMV format image, i.e. we can
        make sense of it. From JH: "The best way to identify images from any of my
        simulators is to look for BEAMLINE=fake in the header."."""

        size, header = FormatSMV.get_smv_header(image_file)

        if header.get("BEAMLINE") == "fake":
            return True
        else:
            return False

    def detectorbase_start(self):
        if not hasattr(self, "detectorbase") or self.detectorbase is None:
            self.detectorbase = SMVImage(self._image_file)
            self.detectorbase.open_file = self.open_file
            self.detectorbase.readHeader()

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header."""

        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector, presuming no one has
        one of these on a two-theta stage. Assert that the beam centre is
        provided in the Mosflm coordinate frame."""

        distance = float(self._header_dictionary["DISTANCE"])
        beam_x = float(self._header_dictionary["BEAM_CENTER_X"])
        beam_y = float(self._header_dictionary["BEAM_CENTER_Y"])
        pixel_size = float(self._header_dictionary["PIXEL_SIZE"])
        image_size = (
            float(self._header_dictionary["SIZE1"]),
            float(self._header_dictionary["SIZE2"]),
        )
        image_pedestal = 1
        try:
            image_pedestal = float(self._header_dictionary["ADC_OFFSET"])
        except (KeyError):
            pass
        overload = 65535 - image_pedestal
        underload = 1 - image_pedestal

        # interpret beam center conventions
        image_height_mm = pixel_size * image_size[1]
        adxv_beam_center = (beam_x, beam_y)
        cctbx_beam_center = (
            adxv_beam_center[0] + pixel_size,
            image_height_mm - adxv_beam_center[1] + pixel_size,
        )

        # Guess whether this is mimicking a Pilatus, if so set detector type so
        # that spot-finding parameters are appropriate
        if pixel_size == 0.172:
            stype = "SENSOR_PAD"
        else:
            stype = "CCD"

        return self._detector_factory.simple(
            stype,
            distance,
            cctbx_beam_center,
            "+x",
            "-y",
            (pixel_size, pixel_size),
            image_size,
            (underload, overload),
            [],
            pedestal=int(self._header_dictionary.get("ADC_OFFSET", 1)),
        )

    def _beam(self):
        """Return a simple model for the beam."""

        wavelength = float(self._header_dictionary["WAVELENGTH"])

        return self._beam_factory.simple(wavelength)

    def _scan(self):
        """Return the scan information for this image."""
        exposure_time = 1
        epoch = None

        # PST, PDT timezones not recognised by default...

        epoch = 0
        try:
            date_str = self._header_dictionary["DATE"]
            date_str = date_str.replace("PST", "").replace("PDT", "")
        except KeyError:
            date_str = ""
        for format_string in ["%a %b %d %H:%M:%S %Y", "%a %b %d %H:%M:%S %Z %Y"]:
            try:
                epoch = calendar.timegm(time.strptime(date_str, format_string))
                break
            except ValueError:
                pass

        # assert(epoch)
        osc_start = float(self._header_dictionary["OSC_START"])
        osc_range = float(self._header_dictionary["OSC_RANGE"])

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        assert len(self.get_detector()) == 1
        panel = self.get_detector()[0]
        image_size = panel.get_image_size()
        return self._get_endianic_raw_data(size=image_size)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVJHSim.understand(arg))
