"""
Helper code to assist with reading TIFF file headers, which are by their
nature binary so we need to mess with things like byte swapping.
"""

import struct
import sys

from dxtbx.format.Format import Format

LITTLE_ENDIAN = 1234
BIG_ENDIAN = 4321


def _tiff_byte_order(four_bytes):
    """Determine the byte order for the file from the magic numbers at the
    very start of the file."""

    if b"II" == four_bytes[:2]:
        if struct.unpack("<H", four_bytes[2:])[0] == 42:
            return LITTLE_ENDIAN
    elif b"MM" == four_bytes[:2]:
        if struct.unpack(">H", four_bytes[2:])[0] == 42:
            return BIG_ENDIAN
    return False


def read_basic_tiff_header(filename):
    """Read the TIFF header (assuming for the moment a 4k header...) and
    return ... something."""

    # things we hope to learn from the vanilla TIFF header

    image_width = None
    image_height = None
    image_depth = None
    header_size = None

    with Format.open_file(filename, "rb") as fh:
        tiff_header = fh.read(1024)
    byte_order = _tiff_byte_order(tiff_header[:4])
    if not byte_order:
        return False

    if byte_order == LITTLE_ENDIAN:
        _I = "<I"
        _H = "<H"
    else:
        _I = ">I"
        _H = ">H"

    offset = struct.unpack(_I, tiff_header[4:8])[0]

    ntags = struct.unpack(_H, tiff_header[offset : offset + 2])[0]
    start = offset + 2

    for j in range(ntags):
        type_desc = struct.unpack(_H, tiff_header[start : start + 2])[0]
        start += 2
        type_type = struct.unpack(_H, tiff_header[start : start + 2])[0]
        start += 2
        # type_size = struct.unpack(_I, tiff_header[start : start + 4])[0]
        start += 4
        if type_type == 4:
            type_offset_or_value = struct.unpack(_I, tiff_header[start : start + 4])[0]
            start += 4
        elif type_type == 3:
            type_offset_or_value = struct.unpack(_H, tiff_header[start : start + 2])[0]
            start += 4

        if type_desc == 256:
            image_width = type_offset_or_value
        elif type_desc == 257:
            image_height = type_offset_or_value
        elif type_desc == 258:
            image_depth = type_offset_or_value
        elif type_desc == 273:
            header_size = type_offset_or_value

    return image_width, image_height, image_depth, header_size, byte_order


def _read_tiff_image_description(tiff_header, byte_order):
    """Search the TIFF header for an image description."""

    # OK then let's get started - and let's assume that the size is > 1 kb

    if byte_order == LITTLE_ENDIAN:
        _I = "<I"
        _H = "<H"
    else:
        _I = ">I"
        _H = ">H"

    offset = struct.unpack(_I, tiff_header[4:8])[0]
    ntags = struct.unpack(_H, tiff_header[offset : offset + 2])[0]
    start = offset + 2

    header_text = None

    for j in range(ntags):
        type_desc = struct.unpack(_H, tiff_header[start : start + 2])[0]
        start += 2
        type_type = struct.unpack(_H, tiff_header[start : start + 2])[0]
        start += 2
        type_size = struct.unpack(_I, tiff_header[start : start + 4])[0]
        start += 4
        if type_type == 4:
            type_offset_or_value = struct.unpack(_I, tiff_header[start : start + 4])[0]
            start += 4
        elif type_type == 3:
            type_offset_or_value = struct.unpack(_H, tiff_header[start : start + 2])[0]
            start += 4
        elif type_type == 2:
            type_offset_or_value = struct.unpack(_I, tiff_header[start : start + 4])[0]
            start += 4

        if type_desc == 270:
            start = type_offset_or_value
            end = type_offset_or_value + type_size
            header_text = tiff_header[start:end].strip()

    return header_text


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        width, height, depth, header, order = read_basic_tiff_header(arg)
        print("(%d x %d) @ %d + %d" % (width, height, depth, header))
        tiff_header = Format.open_file(arg, "rb").read(header)
        text = _read_tiff_image_description(tiff_header, order)
        if text:
            print(text)
        else:
            print("No text found")
