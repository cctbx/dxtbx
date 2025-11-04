"""Format class to recognise images from a ThermoFisher Ceta D detector that
have been converted to an ADSC-like SMV with useful metadata. We
want to override the beam model to produce an unpolarised beam and to set the
detector gain to something sensible."""

from __future__ import annotations

import os

from dxtbx.format.FormatSMVADSC import FormatSMVADSC
from dxtbx.model.beam import Probe


# We have to inherit from FormatSMVADSC rather than FormatSMV, so that we can
# activate this class with an enivronment variable in cases where the
# BEAMLINE or DETECTOR_SN are not set to anything specific (for example files
# that lie and set BEAMLINE=ALS831).
class FormatSMVADSCCetaD(FormatSMVADSC):
    @staticmethod
    def understand(image_file):
        # Allow this class to override FormatSMVADSC with an environment variable
        if "FORCE_SMV_AS_CETAD" in os.environ:
            return True

        # Otherwise recognise specific instruments from the header
        size, header = FormatSMVADSC.get_smv_header(image_file)
        check = header.get("BEAMLINE", "") + header.get("DETECTOR_SN", "")
        if "CETAD" in check:
            return True

        return False

    def _detector(self):
        distance = float(self._header_dictionary["DISTANCE"])
        beam_x = float(self._header_dictionary["BEAM_CENTER_X"])
        beam_y = float(self._header_dictionary["BEAM_CENTER_Y"])
        pixel_size = float(self._header_dictionary["PIXEL_SIZE"])
        image_size = (
            float(self._header_dictionary["SIZE1"]),
            float(self._header_dictionary["SIZE2"]),
        )
        # Ceta has gain of > 26 and saturates at about 8000.0 for binning=1
        # according to Thermo Fisher
        binning = {"1x1": 1, "2x2": 2}.get(self._header_dictionary.get("BIN"), 1)
        gain = float(self._header_dictionary.get("GAIN", 26.0))
        saturation = 8000 * binning**2
        trusted_range = (-1000, saturation)
        pedestal = float(self._header_dictionary.get("IMAGE_PEDESTAL", 0))

        return self._detector_factory.simple(
            "PAD",
            distance,
            (beam_y, beam_x),
            "+x",
            "-y",
            (pixel_size, pixel_size),
            image_size,
            trusted_range,
            [],
            gain=gain,
            pedestal=pedestal,
        )

    def _beam(self):
        """Return a simple model for an unpolarised beam."""

        wavelength = float(self._header_dictionary["WAVELENGTH"])
        return self._beam_factory.make_polarized_beam(
            sample_to_source=(0.0, 0.0, 1.0),
            wavelength=wavelength,
            polarization=(0, 1, 0),
            polarization_fraction=0.5,
            probe=Probe.electron,
        )
