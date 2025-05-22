"""
An implementation of the SMV image reader for Rigaku HyPix-6000 images.
Inherits from FormatSMVRigaku.
"""

from __future__ import annotations

import sys

from boost_adaptbx.boost.python import streambuf
from iotbx.detectors.dtrek import DTREKImage
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.ext import is_big_endian, read_uint32, read_uint32_bs
from dxtbx.format.FormatSMVRigaku import FormatSMVRigaku


class FormatSMVRigakuHypix6000(FormatSMVRigaku):
    """A class for reading SMV format Rigaku HyPix-6000 images, and correctly
    constructing a model for the experiment from this."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like a Rigaku HyPix-6000 SMV format image,
        i.e. we can make sense of it. Essentially that will be if it contains
        all of the keys we are looking for."""

        size, header = FormatSMVRigaku.get_smv_header(image_file)

        wanted_header_items = [
            "DETECTOR_NUMBER",
            "DETECTOR_NAMES",
            "CRYSTAL_GONIO_NUM_VALUES",
            "CRYSTAL_GONIO_NAMES",
            "CRYSTAL_GONIO_UNITS",
            "CRYSTAL_GONIO_VALUES",
            "ROTATION",
            "ROTATION_AXIS_NAME",
            "ROTATION_VECTOR",
            "SOURCE_VECTORS",
            "SOURCE_WAVELENGTH",
            "SOURCE_POLARZ",
            "DIM",
            "SIZE1",
            "SIZE2",
        ]

        if any(item not in header for item in wanted_header_items):
            return False

        detector_prefix = header["DETECTOR_NAMES"].split()[0].strip()

        more_wanted_header_items = [
            "DETECTOR_DIMENSIONS",
            "DETECTOR_SIZE",
            "DETECTOR_VECTORS",
            "GONIO_NAMES",
            "GONIO_UNITS",
            "GONIO_VALUES",
            "GONIO_VECTORS",
            "SPATIAL_BEAM_POSITION",
            "SPATIAL_DISTORTION_VECTORS",
        ]

        if any(
            f"{detector_prefix}{item}" not in header
            for item in more_wanted_header_items
        ):
            return False

        descriptive_items = ["DETECTOR_IDENTIFICATION", "DETECTOR_DESCRIPTION"]

        for header_item in descriptive_items:
            test = f"{detector_prefix}{header_item}"
            if test in header and ("hypix-6000" in header[test].lower()):
                return True

        return False

    def detectorbase_start(self):
        self.detectorbase = DTREKImage(self._image_file)
        self.detectorbase.readHeader()

    def _goniometer(self):
        """Initialize the structure for the goniometer - this will need to
        correctly compose the axes given in the image header. In this case
        this is made rather straightforward as the image header has the
        calculated rotation axis stored in it. We could work from the
        rest of the header and construct a goniometer model."""

        axis = tuple(map(float, self._header_dictionary["ROTATION_VECTOR"].split()))

        return self._goniometer_factory.known_axis(axis)

    def _detector(self):
        """Return a model for the detector, allowing for two-theta offsets
        and the detector position. This will be rather more complex..."""

        detector_name = self._header_dictionary["DETECTOR_NAMES"].split()[0].strip()

        detector_axes = self.get_detector_axes(detector_name)
        detector_x = matrix.col(tuple(detector_axes[:3]))
        detector_y = matrix.col(tuple(detector_axes[3:]))

        # Now map these to real axes
        distortion = self.get_distortion(detector_name)
        detector_fast = detector_x * distortion[0] + detector_y * distortion[1]
        detector_slow = detector_x * distortion[2] + detector_y * distortion[3]

        beam_pixels = self.get_beam_pixels(detector_name)
        pixel_size = self.get_pixel_size(detector_name)
        image_size = self.get_image_size(detector_name)

        detector_origin = -(
            beam_pixels[0] * pixel_size[0] * detector_fast
            + beam_pixels[1] * pixel_size[1] * detector_slow
        )

        gonio_axes = self.get_gonio_axes(detector_name)
        gonio_values = self.get_gonio_values(detector_name)
        gonio_units = self._header_dictionary["%sGONIO_UNITS" % detector_name].split()
        gonio_num_axes = int(
            self._header_dictionary["%sGONIO_NUM_VALUES" % detector_name]
        )

        rotations = []
        translations = []

        for j, unit in enumerate(gonio_units):
            axis = matrix.col(gonio_axes[3 * j : 3 * (j + 1)])
            if unit == "deg":
                rotations.append(
                    axis.axis_and_angle_as_r3_rotation_matrix(gonio_values[j], deg=True)
                )
                translations.append(matrix.col((0.0, 0.0, 0.0)))
            elif unit == "mm":
                rotations.append(
                    matrix.sqr((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
                )
                translations.append(gonio_values[j] * axis)
            else:
                raise RuntimeError("unknown axis unit %s" % unit)

        rotations.reverse()
        translations.reverse()

        for j in range(gonio_num_axes):
            detector_fast = rotations[j] * detector_fast
            detector_slow = rotations[j] * detector_slow
            detector_origin = rotations[j] * detector_origin
            detector_origin = translations[j] + detector_origin

        max_trusted_value = int(float(self._header_dictionary["SATURATED_VALUE"]))
        min_trusted_value = 0

        return self._detector_factory.complex(
            "CCD",
            detector_origin.elems,
            detector_fast.elems,
            detector_slow.elems,
            pixel_size,
            image_size,
            (min_trusted_value, max_trusted_value),
        )

    def _beam(self):
        """Return a simple model for the beam."""

        beam_direction = self.get_beam_direction()
        p_fraction, p_plane = self.get_beam_polarization()

        wavelength = float(self._header_dictionary["SCAN_WAVELENGTH"])

        return self._beam_factory.complex(
            beam_direction, p_fraction, p_plane, wavelength
        )

    def _scan(self):
        """Return the scan information for this image."""

        rotation = self.get_rotation()
        epoch = 0

        exposure_time = rotation[3]
        osc_start = rotation[0]
        osc_range = rotation[2]

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        assert len(self.get_detector()) == 1
        panel = self.get_detector()[0]
        size = panel.get_image_size()
        big_endian = self._header_dictionary["BYTE_ORDER"] == "big_endian"

        with self.open_file(self._image_file, "rb") as fh:
            fh.seek(self._header_size)

            if big_endian == is_big_endian():
                raw_data = read_uint32(streambuf(fh), int(size[0] * size[1]))
            else:
                raw_data = read_uint32_bs(streambuf(fh), int(size[0] * size[1]))

        # note that x and y are reversed here
        raw_data.reshape(flex.grid(size[1], size[0]))
        return raw_data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVRigakuHypix6000.understand(arg))
