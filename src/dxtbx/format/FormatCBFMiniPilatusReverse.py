"""
An implementation of the CBF image reader for Pilatus images for beamlines
with a horizontal, but reversed rotation axis.
"""

from __future__ import annotations

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus


class FormatCBFMiniPilatusReverse(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for beamlines with
    a horizontal reversed rotation axis."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "# Detector" in record and "PILATUS" in record:
                if "S/N XX-XXX" in record or "S/N 60-0123" in record:
                    # S/N 60-0123: SSRF BL18U1
                    # S/N XX-XXX: SSRF BL19U1
                    return True

        return False

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header, though for miniCBF
        there are limited options for this."""

        return self._goniometer_factory.single_axis_reverse()
