from __future__ import annotations

import dxtbx.flumpy as flumpy
from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.model import Detector


class FormatCBFMiniPilatusSPring8_1MSN163(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 125, normally
    at Spring8 BL41XU"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 10-0163" in header
            ):
                return True

        return False

    def __init__(self, image_file, **kwargs):
        self._inverted = False
        super().__init__(image_file, **kwargs)

    def _detector(self) -> Detector:
        # Do inversion of header dictionary here... is intialised in
        # _start so we have to intercept afterwards
        if not self._inverted:
            self._inverted = True
            (
                self._cif_header_dictionary["X-Binary-Size-Fastest-Dimension"],
                self._cif_header_dictionary["X-Binary-Size-Second-Dimension"],
            ) = (
                self._cif_header_dictionary["X-Binary-Size-Second-Dimension"],
                self._cif_header_dictionary["X-Binary-Size-Fastest-Dimension"],
            )
        return super()._detector()

    def get_raw_data(self):
        data = super().get_raw_data()
        # Transpose the data
        ndata = flumpy.from_numpy(flumpy.to_numpy(data).transpose().copy())
        return ndata
