from __future__ import absolute_import, division, print_function

import pytest

from dxtbx.model import Spectrum
from scitbx.array_family import flex


def test_spectrum():
    mean_ev_1 = 12398.4187  # 1.0Å
    wavelengths = 2 / 3, 1, 3 / 2

    for wavelength in wavelengths:
        mean_ev = mean_ev_1 / wavelength
        spectrum = Spectrum(
            flex.double([mean_ev - 50, mean_ev, mean_ev + 50]),
            flex.double([0.5, 1.0, 0.5]),
        )
        assert spectrum.get_weighted_energy_eV() == pytest.approx(mean_ev)
        assert spectrum.get_weighted_wavelength() == pytest.approx(wavelength)

        spectrum = Spectrum.from_dict(spectrum.to_dict())
        assert spectrum.get_weighted_energy_eV() == pytest.approx(mean_ev)
        assert spectrum.get_weighted_wavelength() == pytest.approx(wavelength)
