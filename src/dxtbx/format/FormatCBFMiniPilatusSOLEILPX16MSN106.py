"""Set up for Soleil PX1, with full kappa goniometer"""


import sys

from scitbx.array_family import flex

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus


class FormatCBFMiniPilatusSOLEILPX16MSN106(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 106 @
    Synchrotrol Soleil PX1."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# Detector" in record
                and "PILATUS" in record
                and "S/N 60-0106, Soleil" in record
            ):
                return True

        return False

    def _goniometer(self):
        """Construct a goniometer from the records in the mini CBF header."""

        if (
            "Alpha" in self._cif_header_dictionary
            and "Kappa" in self._cif_header_dictionary
        ):
            # Kappa
            alpha = float(self._cif_header_dictionary["Alpha"].split()[0])
            omega = float(self._cif_header_dictionary["Chi"].split()[0])
            kappa = float(self._cif_header_dictionary["Kappa"].split()[0])
            phi = float(self._cif_header_dictionary["Phi"].split()[0])

            axis = self._cif_header_dictionary["Oscillation_axis"]

            scanaxis = {"OMEGA": "Omega", "PHI": "Phi"}

            assert axis in scanaxis

            # this is the direction the arm points in at datum
            direction = "+z"

            return self._goniometer_factory.make_kappa_goniometer(
                alpha, omega, kappa, phi, direction, scanaxis[axis]
            )

        else:
            # Smargon
            phi = float(self._cif_header_dictionary["Phi"].split()[0])
            chi = float(self._cif_header_dictionary["Chi"].split()[0])
            omega = float(self._cif_header_dictionary["Omega"].split()[0])

            names = flex.std_string(("PHI", "CHI", "OMEGA"))
            axes = flex.vec3_double(((1, 0, 0), (0, 0, -1), (1, 0, 0)))
            angles = flex.double((phi, chi, omega))

            axis = self._cif_header_dictionary["Oscillation_axis"].upper()
            assert axis in names, axis
            scan_axis = flex.first_index(names, axis)

            return self._goniometer_factory.make_multi_axis_goniometer(
                axes, angles, names, scan_axis
            )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusSOLEILPX16MSN106.understand(arg))
