"""Format class for miniCBF files from an ELDICO ED-1 electron diffractometer,
which uses a DECTRIS QUADRO detector (EIGER technology)"""

from __future__ import annotations

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger
from dxtbx.model import SimplePxMmStrategy
from dxtbx.model.beam import Probe


class FormatCBFMiniEigerQuadroED1(FormatCBFMiniEiger):
    @staticmethod
    def understand(image_file):
        header = FormatCBFMiniEiger.get_cbf_header(image_file)

        for record in header.split("\n"):
            if record.startswith("# Detector") and not record.startswith("# Detector_"):
                record = record.lower()
                if not ("quadro" in record or "eiger" in record):
                    return False

            # ED-1 has fixed energy of 160 keV
            if "# wavelength" in record:
                try:
                    wl = float(record.split()[-2])
                except ValueError:
                    return False
                if round(wl, 3) != 0.029:
                    return False

        # If we got this far, check also the mime header to ensure a 512*512 pixel image
        mime_header = ""
        in_binary_format_section = False
        with FormatCBFMiniEiger.open_file(image_file, "rb") as fh:
            for record in fh:
                record = record.decode()
                if "--CIF-BINARY-FORMAT-SECTION--" in record:
                    in_binary_format_section = True
                elif in_binary_format_section and record[0] == "X":
                    mime_header += record
                if in_binary_format_section and len(record.strip()) == 0:
                    # http://sourceforge.net/apps/trac/cbflib/wiki/ARRAY_DATA%20Category
                    #    In an imgCIF file, the encoded binary data begins after
                    #    the empty line terminating the header.
                    break

        for record in mime_header.split("\n"):
            if "Dimension" in record:
                if int(record.split()[-1]) != 512:
                    return False

        return True

    def _goniometer(self):
        """Vertical axis"""

        return self._goniometer_factory.known_axis((0, 1, 0))

    def _beam(self):
        """Ensure an unpolarised beam"""

        wavelength = float(self._cif_header_dictionary["Wavelength"].split()[0])
        return self._beam_factory.make_polarized_beam(
            sample_to_source=(0.0, 0.0, 1.0),
            wavelength=wavelength,
            polarization=(0, 1, 0),
            polarization_fraction=0.5,
            probe=Probe.electron,
        )

    def _detector(self):
        distance = float(self._cif_header_dictionary["Detector_distance"].split()[0])

        beam_xy = (
            self._cif_header_dictionary["Beam_xy"]
            .replace("(", "")
            .replace(")", "")
            .replace(",", "")
            .split()[:2]
        )

        beam_x, beam_y = map(float, beam_xy)

        pixel_xy = (
            self._cif_header_dictionary["Pixel_size"]
            .replace("m", "")
            .replace("x", "")
            .split()
        )

        pixel_x, pixel_y = map(float, pixel_xy)

        if "Silicon" in self._cif_header_dictionary:
            thickness = (
                float(self._cif_header_dictionary["Silicon"].split()[2]) * 1000.0
            )
            material = "Si"
        else:
            thickness = 0.450
            material = "Si"

        # These are always both 512
        nx = int(self._cif_header_dictionary["X-Binary-Size-Fastest-Dimension"])
        ny = int(self._cif_header_dictionary["X-Binary-Size-Second-Dimension"])

        if "Count_cutoff" in self._cif_header_dictionary:
            overload = int(self._cif_header_dictionary["Count_cutoff"].split()[0])
        else:
            overload = 1048576
        if overload == 0:
            # TODO To be checked
            overload = 1048576

        minimum_trusted_value = 0

        try:
            identifier = self._cif_header_dictionary["Detector"].encode()
        except KeyError:
            identifier = "Unknown Eiger"

        detector = self._detector_factory.simple(
            sensor="PAD",
            distance=distance * 1000.0,
            beam_centre=(beam_x * pixel_x * 1000.0, beam_y * pixel_y * 1000.0),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(1000 * pixel_x, 1000 * pixel_y),
            image_size=(nx, ny),
            trusted_range=(minimum_trusted_value, overload),
            mask=[],
        )

        # Here we set specifics: parallax correction and
        # QE correction are effectively disabled by setting the simple
        # pixel-to-millimetre strategy and a very high mu value. Gain is not
        # precisely known. We think it might be around 1.7, but here will leave
        # as 1.0 and expect the user to override at import.
        for panel in detector:
            panel.set_thickness(thickness)
            panel.set_material(material)
            panel.set_identifier(identifier)
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(1e10)

        return detector
