from __future__ import annotations

import h5py
import pytest

import dxtbx.format.Registry as Registry
from dxtbx.format.Format import Format, abstract
from dxtbx.format.FormatHDF5 import FormatHDF5


@abstract
class AbstractFormatTest(Format):
    schemes = ["test"]

    @classmethod
    def understand(cls):
        return True


class FormatTest(AbstractFormatTest):
    @classmethod
    def understand(cls):
        return True


def test_abstract(tmp_path):
    assert AbstractFormatTest.is_abstract()
    assert not FormatTest.is_abstract()

    # Write an empty h5 file so that it won't match
    with h5py.File(tmp_path / "test.h5", "w"):
        pass
    assert FormatHDF5.is_abstract()
    with pytest.raises(TypeError):
        FormatHDF5(str(tmp_path / "test.h5"))

    assert Registry.get_format_class_for_file(str(tmp_path / "test.h5")).is_abstract()
