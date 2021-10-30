"""Byte offset implementation of fullCBF format, for use with CBF files using
byte offset compression which are _not_ made by dectris."""

import binascii
from io import TextIOWrapper
from typing import IO

from dxtbx.ext import uncompress
from dxtbx.format.FormatCBFFull import FormatCBFFull


def _read_cif_binary_format_section(file: IO[bytes]) -> bytes:
    """
    Read a CIF binary format section

    This is a MIME-variant.
    """
    header = ""
    reader = TextIOWrapper(file, encoding="ascii", errors="surrogateescape")
    for line in reader:
        if not line.strip():
            break
        header = header + line
    return header


class FormatCBFFullByteOffset(FormatCBFFull):
    """An image reading class for full CBF format images"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        header = FormatCBFFull.get_cbf_header(image_file)

        # If this is a pilatus cbf file, then we use a different format
        for record in header.split("\n"):
            if "_array_data.header_convention" in record and "PILATUS" in record:
                return False

        with FormatCBFFull.open_file(image_file, "rb") as file:
            # If we've gotten a header, this ends at the MIME header
            file.seek(len(header))
            mime_divider = b"--CIF-BINARY-FORMAT-SECTION--\r\n"
            if not file.read(len(mime_divider)) == mime_divider:
                return False

            binary_header = _read_cif_binary_format_section(file)

        # If we've got this signature in the binary header, we can read
        if 'conversions="x-CBF_BYTE_OFFSET"' in binary_header:
            return True

        return False

    def _read_cbf_image(self):
        start_tag = binascii.unhexlify("0c1a04d5")

        with self.open_file(self._image_file, "rb") as fh:
            data = fh.read()
        data_offset = data.find(start_tag) + 4
        cbf_header = self._parse_cbf_header(
            data[: data_offset - 4].decode("ascii", "ignore")
        )

        pixel_values = uncompress(
            packed=data[data_offset : data_offset + cbf_header["size"]],
            fast=cbf_header["fast"],
            slow=cbf_header["slow"],
        )

        return pixel_values

    def get_raw_data(self):
        if self._raw_data is None:
            data = self._read_cbf_image()
            self._raw_data = data

        return self._raw_data
