from __future__ import annotations

from dxtbx.format import Registry


def test_reading_refl_failure(dials_data):
    test_data = dials_data("centroid_test_data", pathlib=True)

    # Without dials_regression, none of the dials-data tests check for this "invalid binary data" case
    assert Registry.get_format_class_for_file(test_data / "indexed.refl") is None
    # Check .expt while here
    assert Registry.get_format_class_for_file(test_data / "indexed.expt") is None
