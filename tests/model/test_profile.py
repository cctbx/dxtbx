import pytest

from dxtbx.model import ProfileModelFactory


def test_profile_modelling():
    grs = pytest.importorskip("dials.algorithms.profile_model.gaussian_rs")
    if not hasattr(grs, "__file__"):
        pytest.skip("test requires DIALS")
        # this may pretend to be present in Python 3 without dials actually being there

    profile1 = grs.Model(None, 3, 0.1, 0.2, deg=True)
    dictionary = profile1.to_dict()
    profile2 = ProfileModelFactory.from_dict(dictionary)
    assert profile1.sigma_b() == pytest.approx(profile2.sigma_b(), abs=1e-7)
    assert profile1.sigma_m() == pytest.approx(profile2.sigma_m(), abs=1e-7)
