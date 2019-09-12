"""
Experimental format for TIA .ser files used by some FEI microscopes. See
http://www.er-c.org/cbb/info/TIAformat/
"""

from __future__ import absolute_import, division, print_function

import struct
import sys
from builtins import range

from boost.python import streambuf
from scitbx.array_family import flex

import dxtbx.ext
from dxtbx import IncorrectFormatError
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage


class FormatSER(FormatMultiImage, Format):
    def __init__(self, image_file, **kwargs):

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatMultiImage.__init__(self, **kwargs)
        Format.__init__(self, image_file, **kwargs)

    @staticmethod
    def understand(image_file):
        try:
            with FormatSER.open_file(image_file, "rb") as fh:
                tag = fh.read(14)
        except IOError:
            return False

        # File should be little endian
        if tag[0:2] != b"II":
            return False

        # SeriesID: 0x0197 indicates ES Vision Series Data File
        if struct.unpack("<H", tag[2:4])[0] != 0x0197:
            return False

        # SeriesVersion: 0x0210 or 0x0220
        if struct.unpack("<H", tag[4:6])[0] not in (0x0210, 0x0220):
            return False

        # DataTypeID: 0x4122 if elements are 2D arrays
        if struct.unpack("<I", tag[6:10])[0] != 0x4122:
            return False

        # TagTypeID: 0x4142 or 0x4152
        if struct.unpack("<I", tag[10:])[0] not in (0x4142, 0x4152):
            return False

        return True

    @staticmethod
    def _read_metadata(image_file):
        hd = {}
        with FormatSER.open_file(image_file, "rb") as f:
            f.seek(4)
            version = struct.unpack("<H", f.read(2))[0]

            # version specific details
            if version == 528:
                size = 4
                int_fmt = "I"
            else:
                size = 8
                int_fmt = "Q"

            f.seek(14)
            hd["TotalNumberElements"] = struct.unpack("<I", f.read(4))[0]
            hd["ValidNumberElements"] = struct.unpack("<I", f.read(4))[0]
            hd["OffsetArrayOffset"] = struct.unpack("<" + int_fmt, f.read(size))[0]
            hd["NumberDimensions"] = struct.unpack("<I", f.read(4))[0]

            dimension_array = []
            for i in range(hd["NumberDimensions"]):
                d = {}
                d["DimensionSize"] = struct.unpack("<I", f.read(4))[0]
                d["CalibrationOffset"] = struct.unpack("<d", f.read(8))[0]
                d["CalibrationDelta"] = struct.unpack("<d", f.read(8))[0]
                d["CalibrationElement"] = struct.unpack("<I", f.read(4))[0]
                d["DescriptionLength"] = struct.unpack("<I", f.read(4))[0]
                d["Description"] = f.read(d["DescriptionLength"])
                d["UnitsLength"] = struct.unpack("<I", f.read(4))[0]
                d["Units"] = f.read(d["UnitsLength"])
                dimension_array.append(d)
            hd["DimensionArray"] = dimension_array

            f.seek(hd["OffsetArrayOffset"])
            hd["DataOffsetArray"] = struct.unpack(
                "<" + int_fmt * hd["TotalNumberElements"],
                f.read(size * hd["TotalNumberElements"]),
            )
            hd["TagOffsetArray"] = struct.unpack(
                "<" + int_fmt * hd["TotalNumberElements"],
                f.read(size * hd["TotalNumberElements"]),
            )

            # get expected image size by reading metadata for the first image
            f.seek(hd["DataOffsetArray"][0] + 42)
            hd["ArraySizeX"] = struct.unpack("<I", f.read(4))[0]
            hd["ArraySizeY"] = struct.unpack("<I", f.read(4))[0]

        return hd

    def _start(self):
        """Open the image file, read useful metadata into an internal dictionary
        self._header_dictionary"""

        self._header_dictionary = self._read_metadata(self._image_file)

    def get_num_images(self):
        return self._header_dictionary["ValidNumberElements"]

    # This is still required for dials_regression/test.py
    def get_detectorbase(self):
        pass

    def get_goniometer(self, index=None):
        return Format.get_goniometer(self)

    def get_detector(self, index=None):
        return Format.get_detector(self)

    def get_beam(self, index=None):
        return Format.get_beam(self)

    def get_scan(self, index=None):
        if index is None:
            return Format.get_scan(self)
        else:
            scan = Format.get_scan(self)
            return scan[index]

    def get_image_file(self, index=None):
        return Format.get_image_file(self)

    def get_raw_data(self, index):
        data_offset = self._header_dictionary["DataOffsetArray"][index]
        with FormatSER.open_file(self._image_file, "rb") as f:
            f.seek(data_offset)
            d = {}
            d["CalibrationOffsetX"] = struct.unpack("<d", f.read(8))[0]
            d["CalibrationDeltaX"] = struct.unpack("<d", f.read(8))[0]
            d["CalibrationElementX"] = struct.unpack("<I", f.read(4))[0]
            d["CalibrationOffsetY"] = struct.unpack("<d", f.read(8))[0]
            d["CalibrationDeltaY"] = struct.unpack("<d", f.read(8))[0]
            d["CalibrationElementY"] = struct.unpack("<I", f.read(4))[0]
            d["DataType"] = struct.unpack("<H", f.read(2))[0]

            if d["DataType"] == 6:
                read_pixel = dxtbx.ext.read_int32
            elif d["DataType"] == 5:
                read_pixel = dxtbx.ext.read_int16
            elif d["DataType"] == 3:
                read_pixel = dxtbx.ext.read_uint32
            elif d["DataType"] == 2:
                read_pixel = dxtbx.ext.read_uint16
            elif d["DataType"] == 1:
                read_pixel = dxtbx.ext.read_uint8
            else:
                raise RuntimeError(
                    "Image {} data is of an unsupported type".format(index + 1)
                )

            d["ArraySizeX"] = struct.unpack("<I", f.read(4))[0]
            d["ArraySizeY"] = struct.unpack("<I", f.read(4))[0]
            nelts = d["ArraySizeX"] * d["ArraySizeY"]
            raw_data = read_pixel(streambuf(f), nelts)

        # Check image size is as expected (same as the first image)
        if d["ArraySizeX"] != self._header_dictionary["ArraySizeX"]:
            raise RuntimeError(
                "Image {} has an unexpected array size in X".format(index + 1)
            )
        if d["ArraySizeY"] != self._header_dictionary["ArraySizeY"]:
            raise RuntimeError(
                "Image {} has an unexpected array size in Y".format(index + 1)
            )

        image_size = (d["ArraySizeX"], d["ArraySizeY"])
        raw_data.reshape(flex.grid(image_size[1], image_size[0]))

        return raw_data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSER.understand(arg))
