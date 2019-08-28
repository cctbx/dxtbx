#!/usr/bin/env python
# FormatCBF.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Base implementation of CBF formats - which is just really a place holder
# which will tell you whether something is a CBF file (or no.)

from __future__ import absolute_import, division, print_function

import six
from dxtbx.format.Format import Format
from dxtbx import IncorrectFormatError


class FormatCBF(Format):
    """An image reading class for CBF format images i.e. those from Dectris
    amongst others. This is just a first base class which will be used to
    determine whether this is really a CBF file."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        with FormatCBF.open_file(image_file, "rb") as fh:
            return b"###CBF" == fh.read(6)

    @staticmethod
    def get_cbf_header(image_file):
        """Obtain the text section of the header, which is assumed to be
        everything before --CIF-BINARY-FORMAT-SECTION-- - N.B. for reasons
        of simplicity will read the file in 4k chunks."""

        marker = b"--CIF-BINARY-FORMAT-SECTION--"

        with FormatCBF.open_file(image_file, "rb") as fin:
            header = fin.read(4096)
            marker_index = header.find(marker)
            marker_search_start = len(header) - len(marker) + 1
            while marker_index <= 0:
                add = fin.read(4096)
                if not add:
                    # If the marker is not contained in the file then we return the
                    # entire file. This behaviour is enforced by test involving
                    # dials_regression/image_examples/ADSC_CBF/thaumatin_die_M1S5_1_asc_0041.cbf
                    marker_index = None
                    break
                header += add
                marker_index = header.find(marker, marker_search_start)
                marker_search_start += 4096
            header = header[:marker_index]
            if six.PY2:
                return header
            return header.decode("ascii", "ignore")

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        super(FormatCBF, self).__init__(str(image_file), **kwargs)

    @staticmethod
    def _parse_cbf_header(cbf_header):
        header = {
            "fast": 0,
            "slow": 0,
            "length": 0,
            "byte_offset": False,
            "no_compression": False,
        }
        for record in cbf_header.split("\n"):
            if record.startswith("X-Binary-Size-Fastest-Dimension:"):
                header["fast"] = int(record.split()[-1])
            elif record.startswith("X-Binary-Size-Second-Dimension:"):
                header["slow"] = int(record.split()[-1])
            elif record.startswith("X-Binary-Number-of-Elements:"):
                header["length"] = int(record.split()[-1])
            elif record.startswith("X-Binary-Size:"):
                header["size"] = int(record.split()[-1])
            elif "conversions" in record:
                if "x-CBF_BYTE_OFFSET" in record:
                    header["byte_offset"] = True
                elif "x-CBF_NONE" in record:
                    header["no_compression"] = True

        assert header["length"] == header["fast"] * header["slow"]

        return header

    def _start(self):
        """Open the image file, read the image header, copy it into memory
        for future inspection."""

        Format._start(self)

        self._cif_header = FormatCBF.get_cbf_header(self._image_file)

        self._mime_header = ""

        in_binary_format_section = False

        with FormatCBF.open_file(self._image_file, "rb") as fh:
            for record in fh:
                record = record.decode()
                if "--CIF-BINARY-FORMAT-SECTION--" in record:
                    in_binary_format_section = True
                elif in_binary_format_section and record[0] == "X":
                    self._mime_header += record
                if in_binary_format_section and len(record.strip()) == 0:
                    # http://sourceforge.net/apps/trac/cbflib/wiki/ARRAY_DATA%20Category
                    #    In an imgCIF file, the encoded binary data begins after
                    #    the empty line terminating the header.
                    break


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatCBF.understand(arg))
