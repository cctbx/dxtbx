"""
An implementation of the SMV image reader for ADSC images. Inherits from
FormatSMVADSC, customised for beamline 8.3.1 at the ALS where J. Holton uses
two-theta offsets in the vertical direction, as well as idiosyncratic ways
of recording the beam centre... which work fine for ADXV...
"""

import math

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN


class FormatSMVADSCSN926(FormatSMVADSCSN):
    """A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this, for instrument number 926."""

    @staticmethod
    def understand(image_file):
        """Check to see if this is ADSC SN 926."""

        # check this is detector serial number 926 (or 907)

        size, header = FormatSMVADSCSN.get_smv_header(image_file)

        if int(header["DETECTOR_SN"]) not in [926, 907]:
            return False

        return True

    def _detector(self):
        """Return a model for a simple detector, allowing for the installation
        on on a two-theta stage. Assert that the beam centre is provided in
        the Mosflm coordinate frame."""

        distance = float(self._header_dictionary["DISTANCE"])
        if "DENZO_X_BEAM" in self._header_dictionary:
            beam_x = float(self._header_dictionary["DENZO_X_BEAM"])
            beam_y = float(self._header_dictionary["DENZO_Y_BEAM"])
        else:
            beam_x = float(self._header_dictionary["BEAM_CENTER_X"])
            beam_y = float(self._header_dictionary["BEAM_CENTER_Y"])
        pixel_size = float(self._header_dictionary["PIXEL_SIZE"])
        image_size = (
            float(self._header_dictionary["SIZE1"]),
            float(self._header_dictionary["SIZE2"]),
        )
        if "TWOTHETA" in self._header_dictionary:
            two_theta = float(self._header_dictionary["TWOTHETA"])
        else:
            two_theta = 0

        # now correct for some idiosyncracies...
        # two-theta included in beam centre - so remove this
        beam_y += distance * math.tan(two_theta * math.pi / 180.0)

        return self._detector_factory.two_theta(
            "CCD",
            distance,
            (beam_y, beam_x),
            "+x",
            "-y",
            "+x",
            two_theta,
            (pixel_size, pixel_size),
            image_size,
            self._adsc_trusted_range(),
            [],
            gain=self._adsc_module_gain(),
        )
