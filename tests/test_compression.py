from scitbx.array_family import flex

from dxtbx.ext import compress, uncompress


def test_compress_decompress():
    x, y = 10, 10

    original = [1 for n in range(x * y)]

    original[10] = 44369
    original[11] = 214

    data = flex.int(original)
    compressed = compress(data)
    uncompressed = uncompress(compressed, x, y)

    assert list(data) == list(uncompressed)
