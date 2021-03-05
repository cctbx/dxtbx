from scitbx.array_family import flex

from dxtbx.ext import compress, uncompress


def test_compress_decompress():
    x, y = 10, 10

    data = flex.int(x * y, 1)
    data[10] = 44369
    data[11] = 214
    compressed = compress(data)
    uncompressed = uncompress(compressed, x, y)

    assert list(data) == list(uncompressed)
