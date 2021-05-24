"""
An implementation of the SMV image reader for ADSC images. Inherits from
FormatSMVADSC, customised for example on ALS beamline 8.3.1 from back in the
day which had its own way of recording beam centre.
"""


import time

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN


class FormatSMVADSCSN442(FormatSMVADSCSN):
    """A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this, for instrument number 442."""

    @staticmethod
    def understand(image_file):
        """Check to see if this is ADSC SN 442."""

        # check this is detector serial number 442

        size, header = FormatSMVADSCSN.get_smv_header(image_file)

        return int(header["DETECTOR_SN"]) == 442

    def _detector(self):
        """Return a model for a simple detector, presuming no one has
        one of these on a two-theta stage. Assert that the beam centre is
        provided in the Mosflm coordinate frame."""

        # this detector lived two lives - it moved on May 14, 2007 to ALS 5.0.1
        # where the denzo beam centre was removed

        date_str = self._header_dictionary["DATE"]
        format_string = "%a %b %d %H:%M:%S %Y"
        tm = time.strptime(date_str, format_string)

        als831 = True
        if tm.tm_year > 2006 or (tm.tm_year == 2006 and tm.tm_yday >= 134):
            als831 = False

        if als831:
            beam_x = float(self._header_dictionary["DENZO_X_BEAM"])
            beam_y = float(self._header_dictionary["DENZO_Y_BEAM"])
        else:
            beam_x = float(self._header_dictionary["BEAM_CENTER_X"])
            beam_y = float(self._header_dictionary["BEAM_CENTER_Y"])

        distance = float(self._header_dictionary["DISTANCE"])
        pixel_size = float(self._header_dictionary["PIXEL_SIZE"])
        image_size = (
            float(self._header_dictionary["SIZE1"]),
            float(self._header_dictionary["SIZE2"]),
        )
        trusted_range = self._adsc_trusted_range(pedestal=40)

        return self._detector_factory.simple(
            "CCD",
            distance,
            (beam_y, beam_x),
            "+x",
            "-y",
            (pixel_size, pixel_size),
            image_size,
            trusted_range,
            [],
            gain=self._adsc_module_gain(),
            pedestal=40,
        )
