from __future__ import absolute_import, division, print_function

import pytest

from dxtbx.model import Spectrum
from scitbx.array_family import flex


def test_spectrum():
    mean_ev = 12398.4187  # 1.0Å

    spectrum = Spectrum(
        flex.double([mean_ev - 50, mean_ev, mean_ev + 50]), flex.double([0.5, 1.0, 0.5])
    )
    assert spectrum.get_weighted_wavelength() == pytest.approx(1.0)

    spectrum = Spectrum.from_dict(spectrum.to_dict())
    assert spectrum.get_weighted_wavelength() == pytest.approx(1.0)
