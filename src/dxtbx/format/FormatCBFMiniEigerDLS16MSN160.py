import sys

import libtbx
from scitbx.array_family import flex

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger
from dxtbx.masking import GoniometerMaskerFactory


class FormatCBFMiniEigerDLS16MSN160(FormatCBFMiniEiger):
    @staticmethod
    def understand(image_file):
        """Check to see if this format class can understand the image file.

        Args:
          image_file (str): The file path of the image file to check.

        Returns:
          bool: Returns ``True`` if the image_file is understood by this format class,
          else returns ``False``.

        """

        # this depends on DIALS for the goniometer shadow model; if missing
        # simply return False

        header = FormatCBFMiniEiger.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# detector" in record.lower()
                and "eiger" in record.lower()
                and "S/N 160-0001" in header
            ):
                return True

        return False

    @staticmethod
    def has_dynamic_shadowing(**kwargs):
        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if dynamic_shadowing in (libtbx.Auto, "Auto"):
            return True
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)
        super().__init__(image_file, **kwargs)

    def _goniometer(self):
        """Return a model for a multi-axis goniometer.

        This should be checked against the image header, though for miniCBF
        there are limited options for this.

        Returns:
          dxtbx.model.Goniometer.MultiAxisGoniometer: The goniometer model for
          this detector.

        """

        if "Phi" in self._cif_header_dictionary:
            phi = float(self._cif_header_dictionary["Phi"].split()[0])
        else:
            phi = 0.0

        if "Chi" in self._cif_header_dictionary:
            chi = float(self._cif_header_dictionary["Chi"].split()[0])
        else:
            chi = 0.0

        if "Omega" in self._cif_header_dictionary:
            omega = float(self._cif_header_dictionary["Omega"].split()[0])
        else:
            omega = 0.0

        phi_axis = (1, 0, 0)
        chi_axis = (0, 0, -1)
        omega_axis = (1, 0, 0)
        axes = flex.vec3_double((phi_axis, chi_axis, omega_axis))
        angles = flex.double((phi, chi, omega))
        names = flex.std_string(("GON_PHI", "GON_CHI", "GON_OMEGA"))
        return self._goniometer_factory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    def get_goniometer_shadow_masker(self, goniometer=None):
        if goniometer is None:
            goniometer = self.get_goniometer()

        assert goniometer is not None

        if goniometer.get_names()[1] == "GON_CHI":
            # SmarGon
            return GoniometerMaskerFactory.smargon(goniometer)

        else:
            raise ValueError(
                "Don't understand this goniometer: %s" % list(goniometer.get_names())
            )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatCBFMiniEigerDLS16MSN160.understand(arg))
