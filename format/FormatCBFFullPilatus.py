#!/usr/bin/env python
# FormatCBFFullPilatus.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Pilatus implementation of fullCBF format, for use with Dectris detectors.

from __future__ import absolute_import, division, print_function

import re

from dxtbx.format.FormatCBFFull import FormatCBFFull
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask


class FormatCBFFullPilatus(FormatCBFFull):
    """An image reading class for full CBF format images from Pilatus
    detectors."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        header = FormatCBFFull.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "_array_data.header_convention" in record and "PILATUS" in record:
                return True

        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatCBFFull.__init__(self, image_file, **kwargs)

        self._raw_data = None

    def _start(self):
        """Open the image file as a cbf file handle, and keep this somewhere
        safe."""

        FormatCBFFull._start(self)

    def _beam(self):
        """Return a working beam instance. Override polarization to be 0.999."""

        beam = self._beam_factory.imgCIF_H(self._get_cbf_handle())
        beam.set_polarization_fraction(0.999)
        return beam

    def _detector(self):
        """Return a working detector instance, with added mask regions."""

        detector = self._detector_factory.imgCIF_H(self._get_cbf_handle(), "PAD")

        for f0, f1, s0, s1 in determine_pilatus_mask(detector):
            detector[0].add_mask(f0 - 1, s0 - 1, f1, s1)

        m = re.search(
            r"^#\s*(\S+)\ssensor, thickness\s*([0-9.]+)\s*m\s*$",
            self._cif_header,
            re.MULTILINE,
        )
        if m:
            # header gives thickness in metres, we store mm
            thickness = float(m.group(2)) * 1000
            material = m.group(1)

            if material == "Silicon":
                material = "Si"

            for panel in detector:
                panel.set_thickness(thickness)
                panel.set_material(material)

            try:
                # a header only CBF file will not have a beam object
                beam = self._beam()

            except Exception:
                beam = None

            if beam:
                # attenuation coefficient depends on the beam wavelength
                wavelength = beam.get_wavelength()

                from cctbx.eltbx import attenuation_coefficient
                from dxtbx.model import ParallaxCorrectedPxMmStrategy

                # this will fail for undefined composite materials
                table = attenuation_coefficient.get_table(material)
                # mu_at_angstrom returns cm^-1
                mu = table.mu_at_angstrom(wavelength) / 10.0

                for panel in detector:
                    panel.set_px_mm_strategy(
                        ParallaxCorrectedPxMmStrategy(mu, thickness)
                    )
                    panel.set_mu(mu)

        m = re.search(r"^#\s*Detector:\s+(.*?)\s*$", self._cif_header, re.MULTILINE)
        if m and m.group(1):
            panel.set_identifier(m.group(1).encode())

        size = detector[0].get_image_size()
        if size == (2463, 2527):
            self.vendortype = "Pilatus-6M"
        elif size == (1475, 1679):
            self.vendortype = "Pilatus-2M"
        elif size == (487, 619):
            self.vendortype = "Pilatus-300K"

        return detector

    def _read_cbf_image(self):
        from cbflib_adaptbx import uncompress
        import binascii

        start_tag = binascii.unhexlify("0c1a04d5")

        with self.open_file(self._image_file, "rb") as fh:
            data = fh.read()
        data_offset = data.find(start_tag) + 4
        cbf_header = self._parse_cbf_header(
            data[: data_offset - 4].decode("ascii", "ignore")
        )

        pixel_values = uncompress(
            packed=data[data_offset : data_offset + cbf_header["size"]],
            fast=cbf_header["fast"],
            slow=cbf_header["slow"],
        )

        return pixel_values

    def get_raw_data(self):
        if self._raw_data is None:
            data = self._read_cbf_image()
            self._raw_data = data

        return self._raw_data

    def get_vendortype(self):
        from dxtbx.format.FormatPilatusHelpers import get_vendortype as gv

        return gv(self.get_detector())


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatCBFFullPilatus.understand(arg))
