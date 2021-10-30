"""
An implementation of the CBF image reader for Pilatus images, from the Pilatus
6M SN 100 currently on Diamond I04.
"""


from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus


class FormatCBFMiniPilatusXXX(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN XXX."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N XX-XXX" in header
            ):
                return True

        return False

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header, though for miniCBF
        there are limited options for this."""

        return self._goniometer_factory.single_axis_reverse()
