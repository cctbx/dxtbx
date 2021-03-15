import h5py
import numpy
import pytest

from scitbx.array_family import flex  # noqa: F401; boost python bindings

from dxtbx_format_nexus_ext import (
    dataset_as_flex_double,
    dataset_as_flex_float,
    dataset_as_flex_int,
)


def test_dataset_as_flex_int():
    # Create an in-memory HDF5 dataset
    f = h5py.File("int.h5", "w", driver="core", backing_store=False)

    shape = (20, 20, 20)
    dataset = f.create_dataset(
        "data",
        shape,
        dtype="int",
    )

    # create some random data
    array = numpy.random.rand(*shape) * 100
    original = array.astype("int")

    dataset[:] = original

    selection = slice(0, shape[0], 1), slice(0, shape[1], 1), slice(0, shape[2], 1)
    foo = dataset_as_flex_int(dataset.id.id, selection)

    assert foo.as_numpy_array() == pytest.approx(original)
    f.close()


def test_dataset_as_flex_float():
    # Create an in-memory HDF5 dataset
    f = h5py.File("float.h5", "w", driver="core", backing_store=False)

    shape = (20, 20, 20)
    dataset = f.create_dataset(
        "data",
        shape,
        dtype="float",
    )

    # create some random data
    array = numpy.random.rand(*shape) * 100
    original = array.astype("float")

    dataset[:] = original

    selection = slice(0, shape[0], 1), slice(0, shape[1], 1), slice(0, shape[2], 1)
    foo = dataset_as_flex_float(dataset.id.id, selection)

    assert foo.as_numpy_array() == pytest.approx(original)
    f.close()


def test_dataset_as_flex_double():
    # Create an in-memory HDF5 dataset
    f = h5py.File("double.h5", "w", driver="core", backing_store=False)

    shape = (20, 20, 20)
    dataset = f.create_dataset(
        "data",
        shape,
        dtype="double",
    )

    # create some random data
    array = numpy.random.rand(*shape) * 100
    original = array.astype("double")

    dataset[:] = original

    selection = slice(0, shape[0], 1), slice(0, shape[1], 1), slice(0, shape[2], 1)
    foo = dataset_as_flex_double(dataset.id.id, selection)

    assert foo.as_numpy_array() == pytest.approx(original)
    f.close()
