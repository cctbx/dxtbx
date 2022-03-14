from __future__ import annotations

import h5py
import hdf5plugin
import numpy
import pytest

from scitbx.array_family import flex  # noqa: F401; boost python bindings

from dxtbx_format_nexus_ext import (
    dataset_as_flex_double,
    dataset_as_flex_float,
    dataset_as_flex_int,
)


def uncompressed(file, shape, type_name):
    return file.create_dataset(
        "data",
        shape,
        dtype=type_name,
    )


def gzip(file, shape, type_name):
    return file.create_dataset(
        "data",
        shape,
        compression="gzip",
        dtype=type_name,
    )


def bshuf_lz4(file, shape, type_name):
    return file.create_dataset(
        "data",
        shape,
        **hdf5plugin.Bitshuffle(),
        dtype=type_name,
    )


@pytest.mark.parametrize(
    "creator",
    [
        uncompressed,
        gzip,
        bshuf_lz4,
    ],
)
@pytest.mark.parametrize(
    "type_name,converter",
    [
        ("int", dataset_as_flex_int),
        ("float", dataset_as_flex_float),
        ("double", dataset_as_flex_double),
    ],
)
def test_dataset_as_flex(type_name, creator, converter):
    # Create an in-memory HDF5 dataset with unique name
    f = h5py.File(type_name + ".h5", "w", driver="core", backing_store=False)

    shape = (20, 20, 20)
    dataset = creator(f, shape, type_name)

    # create some random data
    array = numpy.random.rand(*shape) * 100
    original = array.astype(type_name)

    dataset[:] = original

    selection = slice(0, shape[0], 1), slice(0, shape[1], 1), slice(0, shape[2], 1)
    foo = converter(dataset.id.id, selection)

    assert foo.as_numpy_array() == pytest.approx(original)
