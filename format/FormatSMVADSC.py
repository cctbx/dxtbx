"""An implementation of the SMV image reader for ADSC images"""


import calendar
import sys
import time

from iotbx.detectors import SMVImage

from dxtbx.format.FormatSMV import FormatSMV


class FormatSMVADSC(FormatSMV):
    """A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an ADSC SMV format image, i.e. we
        can make sense of it. Essentially that will be if it contains all of
        the keys we are looking for and not some we are not (i.e. that belong
        to a Rigaku Saturn.)"""

        size, header = FormatSMV.get_smv_header(image_file)

        # do not understand JHSim images
        if header.get("BEAMLINE") == "fake":
            return False

        # do not understand Timepix_SU images
        if header.get("BEAMLINE", "").upper() == "TIMEPIX_SU":
            return False

        # this used to include TIME
        wanted_header_items = [
            "BEAM_CENTER_X",
            "BEAM_CENTER_Y",
            "DISTANCE",
            "WAVELENGTH",
            "PIXEL_SIZE",
            "OSC_START",
            "OSC_RANGE",
            "SIZE1",
            "SIZE2",
            "BYTE_ORDER",
        ]

        if any(item not in header for item in wanted_header_items):
            return False

        unwanted_header_items = ["DTREK_DATE_TIME"]

        if any(item in header for item in unwanted_header_items):
            return False

        return True

    def detectorbase_start(self):
        if not hasattr(self, "detectorbase") or self.detectorbase is None:
            self.detectorbase = SMVImage(self._image_file)
            self.detectorbase.open_file = self.open_file
            self.detectorbase.readHeader()

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. Invert the axis if
        the rotation range is negative."""

        if float(self._header_dictionary["OSC_RANGE"]) > 0:
            return self._goniometer_factory.single_axis()
        else:
            return self._goniometer_factory.single_axis_reverse()

    def _adsc_trusted_range(self, pedestal=None):
        """Return a 16 bit trusted range shifted to account for any image
        pedestal that is present"""

        if pedestal is None:
            pedestal = float(self._header_dictionary.get("IMAGE_PEDESTAL", 0))

        overload = 65535 - pedestal
        underload = -1 - pedestal

        return underload, overload

    def _adsc_module_gain(self, model=None):
        """Return an appropriate gain value in ADU per captured X-ray for an
        ADSC CCD module. If the model is None (unknown) then make a guess based
        on header details"""

        # Values are based on a combination of manufacturer datasheets,
        # James Holton's list at http://bl831.als.lbl.gov/xtalsize.html and
        # empirical guesswork based on spot finding results. The accuracy should
        # not be expected to be better than 20% and indeed may be worse than that.
        # Values may refer to gain in ADU per incident photon rather than the
        # more appropriate ADU per captured photon.

        # Allow gain set in the header to override the guesswork here
        if "GAIN" in self._header_dictionary:
            return float(self._header_dictionary["GAIN"])

        # Get the binning level
        bin_lev = str(self._header_dictionary.get("BIN"))

        if model is None:
            # guesswork based on the header details
            npx1 = int(self._header_dictionary["SIZE1"])
            npx2 = int(self._header_dictionary["SIZE2"])
            px_sz = float(self._header_dictionary["PIXEL_SIZE"])
            if npx1 == npx2 == 2304 and px_sz == 0.0816:
                model = "Q4"  # or Q4R
            elif npx1 == npx2 == 4096 and px_sz == 0.0512:
                model = "Q210"  # or Q210R bin none or 1x1
            elif npx1 == npx2 == 2048 and px_sz == 0.1024:
                model = "Q210"  # or Q210R bin 2x2
                bin_lev = "2x2"
            elif npx1 == npx2 == 6144 and px_sz == 0.0512:
                model = "Q315"  # or Q315R bin none or 1x1
            elif npx1 == npx2 == 3072 and px_sz == 0.1024:
                model = "Q315"  # or Q315R bin 2x2
                bin_lev = "2x2"
            elif npx1 == npx2 == 4168 and px_sz == 0.0648:
                model = "Q270"  # or Q270R bin none or 1x1
            elif npx1 == npx2 == 2084 and px_sz == 0.0324:
                model = "Q270"  # or Q270R bin 2x2
                bin_lev = "2x2"
            else:
                # Unidentified model
                return 1.0

        model = model.upper()
        if model.startswith("Q270"):
            # Don't know the difference between binning modes in this case
            return 2.8

        if model.startswith("Q4"):
            # https://web.archive.org/web/20051224172252/http://www.adsc-xray.com:80/prod_compare.html
            # Don't know the difference between binning modes in this case
            return 1.25

        # Get binning mode HW/SW
        bin_type = self._header_dictionary.get("BIN_TYPE")

        if bin_lev.upper() == "NONE":
            bin_type = None
        elif bin_type != "HW":
            bin_type = "SW"

        if bin_type == "SW":
            # return 0.6 This does not look believable. Default to 2.4 instead
            return 2.4

        if bin_type != "HW":  # assume unbinned
            return 2.4
        elif model.endswith("R"):  # e.g. Q315r, HW binning
            return 1.8
        else:  # e.g. Q315, HW binning
            return 2.4

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

        return self._detector_factory.simple(
            "CCD",
            distance,
            (beam_y, beam_x),
            "+x",
            "-y",
            (pixel_size, pixel_size),
            image_size,
            self._adsc_trusted_range(),
            [],
            gain=self._adsc_module_gain(),
            pedestal=float(self._header_dictionary.get("IMAGE_PEDESTAL", 0)),
        )

    def _beam(self):
        """Return a simple model for the beam."""

        wavelength = float(self._header_dictionary["WAVELENGTH"])

        return self._beam_factory.simple(wavelength)

    def _scan(self):
        """Return the scan information for this image."""
        exposure_time = float(self._header_dictionary["TIME"])
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
        osc_range = abs(float(self._header_dictionary["OSC_RANGE"]))

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
        print(FormatSMVADSC.understand(arg))
