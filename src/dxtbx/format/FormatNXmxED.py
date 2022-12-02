from __future__ import annotations

import h5py

from scitbx.array_family import flex

import dxtbx.nexus
from dxtbx.format import nexus
from dxtbx.format.FormatNXmx import FormatNXmx
from dxtbx.masking import mask_untrusted_circle, mask_untrusted_polygon
from dxtbx.model import SimplePxMmStrategy


class FormatNXmxED(FormatNXmx):
    """Format class for 3DED/MicroED images that have been converted to NXmx
    by nexgen"""

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            if FormatNXmxED.get_source_probe(handle).lower() != "electron":
                return False
        return True

    @staticmethod
    def get_source_probe(handle) -> str:
        try:
            probe = handle["/entry/source/probe"][()]
            probe = nexus.h5str(probe)
        except KeyError:
            probe = ""
        return probe

    def _detector(self):
        """Overridden to change px-to-mm and QE corrections from their X-ray
        defaults"""

        # DECTRIS recommend a paper by Fernandez-Perez et al 2021
        # (https://doi.org/10.1088/1748-0221/16/10/P10034) that provides some
        # characterisation of the detector response to fast electrons. Fig 1
        # shows there is essentially no transmission through the sensor layer.
        # Most charge is accumulated near the surface of the sensor, but due
        # to the random-walk effect of multiple interactions as the electron
        # slows down, there is significant spreading of the signal beyond the
        # incident pixel. The parallax-corrected px-to-mm function used for
        # X-rays is not correct for this situation, so we will set the simple
        # strategy instead.
        #
        # A second issue is that of the QE correction. As there is no
        # transmission of electrons through the sensor, the simple model using
        # a linear absorption coefficient, µ, is also not appropriate here.
        # There is loss of signal due to other effects, such as backscattering,
        # but we do not have a model for this. As a result, we shall set µ to be
        # very high to effectively nullify the correction.

        detector = super()._detector()
        for panel in detector:
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(1e10)

        return detector

    def _beam(self, index=None):
        """Ensure the beam is unpolarised"""

        beam = super()._beam()
        beam.set_polarization_fraction(0.5)

        return beam


class FormatNXmxEDeBIC(FormatNXmxED):

    """Set a mask specific for the SINGLA installed at eBIC"""

    @staticmethod
    def understand(image_file):
        with h5py.File(image_file) as handle:
            name = dxtbx.nexus.nxmx.h5str(FormatNXmxEDeBIC.get_instrument_name(handle))
            if name and "ebic" in name.lower():
                return True
        return False

    def _start(self):
        super()._start()

        self._extend_static_mask()

    def _extend_static_mask(self):
        """Extend the existing static mask to include the flange shadow."""

        pixel_mask = self._static_mask[0]

        # It is easier to calculate the inverted mask for the shadow of the vacuum flange
        inverted_flange_mask = flex.bool(flex.grid(pixel_mask.all()), True)
        for coord in [(201, 210), (202, 846), (836, 208), (838, 845)]:
            mask_untrusted_circle(inverted_flange_mask, coord[0], coord[1], 194)
        vertices = [
            (7, 210),
            (201, 16),
            (836, 14),
            (1028, 174),
            (1028, 886),
            (838, 1039),
            (202, 1040),
            (8, 836),
        ]
        polygon = flex.vec2_double(vertices)
        mask_untrusted_polygon(inverted_flange_mask, polygon)

        # Combine masks taking care of the inversion
        self._static_mask = ((pixel_mask & ~inverted_flange_mask),)
