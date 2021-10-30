"""Byte offset implementation of fullCBF format, for use with CBF files using
byte offset compression which are _not_ made by dectris."""

from __future__ import absolute_import, division, print_function

import binascii

from dxtbx.ext import uncompress
from dxtbx.format.FormatCBFFull import FormatCBFFull


class FormatCBFFullByteOffset(FormatCBFFull):
    """An image reading class for full CBF format images."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        header = FormatCBFFull.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "_array_data.header_convention" in record and "PILATUS" in record:
                return False

        # this is ugly but I don't know a better way to seek for this text
        if b'conversions="x-CBF_BYTE_OFFSET"' in FormatCBFFull.open_file(
            image_file, "rb"
        ).read(10000):
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
