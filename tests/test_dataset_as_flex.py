import h5py
import numpy
import pytest

from scitbx.array_family import flex  # noqa: F401; boost python bindings

from dxtbx_format_nexus_ext import (
    dataset_as_flex_double,
    dataset_as_flex_float,
    dataset_as_flex_int,
)

try:
    import bitshuffle.h5
except ImportError:
    bitshuffle = None


def uncompressed(file, shape, type_name):
    return file.create_dataset(
        "data",
        shape,
        dtype=type_name,
    )


def bshuf_lz4(file, shape, type_name):
    block_size = 0
    return file.create_dataset(
        "data",
        shape,
        compression=bitshuffle.h5.H5FILTER,
        compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
        dtype=type_name,
    )


@pytest.mark.parametrize(
    "creator",
    [
        uncompressed,
        pytest.param(
            bshuf_lz4,
            marks=pytest.mark.skipif(
                bitshuffle is None, reason="bitshuffle module not available"
            ),
        ),
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
