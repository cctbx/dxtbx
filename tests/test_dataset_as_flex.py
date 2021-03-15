import h5py
import numpy
import pytest

from scitbx.array_family import flex  # noqa: F401; boost python bindings

from dxtbx_format_nexus_ext import (
    dataset_as_flex_double,
    dataset_as_flex_float,
    dataset_as_flex_int,
)


@pytest.mark.parametrize(
    "type_name,converter",
    [
        ("int", dataset_as_flex_int),
        ("float", dataset_as_flex_float),
        ("double", dataset_as_flex_double),
    ],
)
def test_dataset_as_flex(type_name, converter):
    # Create an in-memory HDF5 dataset with unique name
    f = h5py.File(type_name + ".h5", "w", driver="core", backing_store=False)

    shape = (20, 20, 20)
    dataset = f.create_dataset(
        "data",
        shape,
        dtype=type_name,
    )

    # create some random data
    array = numpy.random.rand(*shape) * 100
    original = array.astype(type_name)

    dataset[:] = original

    selection = slice(0, shape[0], 1), slice(0, shape[1], 1), slice(0, shape[2], 1)
    foo = converter(dataset.id.id, selection)

    assert foo.as_numpy_array() == pytest.approx(original)
