import pytest
import six.moves.cPickle as pickle

from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.model import (
    Crystal,
    CrystalFactory,
    MosaicCrystalKabsch2010,
    MosaicCrystalSauter2014,
)


@pytest.fixture()
def example_crystal():
    real_space_a = matrix.col((35.2402102454, -7.60002142787, 22.080026774))
    real_space_b = matrix.col((22.659572494, 1.47163505925, -35.6586361881))
    real_space_c = matrix.col((5.29417246554, 38.9981792999, 4.97368666613))
    return {
        "real_space_a": real_space_a,
        "real_space_b": real_space_b,
        "real_space_c": real_space_c,
        "space_group_symbol": "P 1 2/m 1",
    }


def test_crystal(example_crystal):
    c1 = Crystal(**example_crystal)

    d = c1.to_dict()
    c2 = CrystalFactory.from_dict(d)
    for direction in ("real_space_a", "real_space_b", "real_space_c"):
        assert abs(matrix.col(d[direction]) - example_crystal[direction]) <= 1e-7
    assert d["space_group_hall_symbol"] == "-P 2y"
    assert c1 == c2


def test_mosaic_crystal(example_crystal):
    c1 = MosaicCrystalKabsch2010(**example_crystal)
    c1.set_mosaicity(0.1)

    d = c1.to_dict()
    c2 = CrystalFactory.from_dict(d)
    for direction in ("real_space_a", "real_space_b", "real_space_c"):
        assert abs(matrix.col(d[direction]) - example_crystal[direction]) <= 1e-7
    assert d["space_group_hall_symbol"] == "-P 2y"
    assert d["mosaicity"] == 0.1
    assert c1 == c2


def test_crystal_with_scan_points(example_crystal):
    c1 = Crystal(**example_crystal)

    A = c1.get_A()
    c1.set_A_at_scan_points([A for i in range(5)])

    # Set the B covariance. The values are nonsense, just ensure they are
    # all different
    cov_B = flex.double(range(9 * 9)) * 1e-5
    c1.set_B_covariance(cov_B)
    cov_B.reshape(flex.grid(1, 9, 9))
    cov_B_array = flex.double(flex.grid(5, 9, 9))
    for i in range(5):
        cov_B_array[i : (i + 1), :, :] = cov_B
    c1.set_B_covariance_at_scan_points(cov_B_array)
    cov_B = c1.get_B_covariance()

    d = c1.to_dict()
    c2 = CrystalFactory.from_dict(d)
    eps = 1e-9
    for Acomp in d["A_at_scan_points"]:
        for e1, e2 in zip(A, Acomp):
            assert abs(e1 - e2) <= eps
    for covBcomp in d["B_covariance_at_scan_points"]:
        for e1, e2 in zip(cov_B, covBcomp):
            assert abs(e1 - e2) <= eps

    assert c1 == c2


@pytest.mark.parametrize(
    "crystal_class", [Crystal, MosaicCrystalKabsch2010, MosaicCrystalSauter2014]
)
def test_crystal_with_recalculated_cell(crystal_class, example_crystal):
    c1 = crystal_class(**example_crystal)
    uc = c1.get_unit_cell()
    c1.set_recalculated_unit_cell(uc)
    c1.set_recalculated_cell_parameter_sd((0.1,) * 6)
    c1.set_recalculated_cell_volume_sd(0.001)

    d = c1.to_dict()
    c2 = CrystalFactory.from_dict(d)
    c3 = pickle.loads(pickle.dumps(c1))

    for c in (c2, c3):
        assert c.get_recalculated_unit_cell() is not None
        assert c1.get_recalculated_unit_cell().is_similar_to(
            c.get_recalculated_unit_cell()
        )
        assert c1 == c
        assert c.get_recalculated_cell_parameter_sd() == (0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
        assert c.get_recalculated_cell_volume_sd() == 0.001
