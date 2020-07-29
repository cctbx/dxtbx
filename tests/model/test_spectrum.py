from __future__ import absolute_import, division, print_function

import math
import random

import pytest

from scitbx.array_family import flex

from dxtbx.model import Spectrum


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


def test_spectrum_bandwidth():
    e, c = __generate_spectrum()

    spectrum = Spectrum(e, c)

    emin, emax = spectrum.get_emin_ev(), spectrum.get_emax_ev()

    assert emin == pytest.approx(11840, abs=2)
    assert emax == pytest.approx(11876, abs=2)


def __generate_spectrum():
    energies = flex.double()
    counts = flex.double()

    peaks, heights, widths = (
        (11860, 11840, 11853, 11860, 11865, 11870),
        (10000, 5000, 5000, 5000, 5000),
        (10, 1, 1, 1, 1),
    )

    for dev in range(118000, 120000):
        ev = 0.1 * dev
        c = random.random() * 100 - 50
        for p, h, w in zip(peaks, heights, widths):
            c += h * math.exp(-(((ev - p) / w) ** 2))
        energies.append(ev)
        counts.append(c)

    return energies, counts
