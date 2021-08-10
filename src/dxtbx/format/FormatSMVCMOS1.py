"""An implementation of the SMV image reader for CMOS1 images, from ALS 4.2.2"""


import calendar
import sys
import time

from boost_adaptbx.boost.python import streambuf
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.ext import read_uint16
from dxtbx.format.FormatSMV import FormatSMV


class FormatSMVCMOS1(FormatSMV):
    """A class for reading SMV format CMOS1 images. 'Abstract' class."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like a CMOS1 d*TREK SMV format image,
        i.e. we can make sense of it. Essentially that will be if it contains
        all of the keys we are looking for."""

        size, header = FormatSMV.get_smv_header(image_file)

        wanted_header_items = [
            "DETECTOR_NUMBER",
            "DETECTOR_NAMES",
            "BYTE_ORDER",
            "DIM",
            "SIZE1",
            "SIZE2",
            "Data_type",
        ]

        if any(item not in header for item in wanted_header_items):
            return False

        detector_prefixes = header["DETECTOR_NAMES"].split()

        if len(detector_prefixes) != 1:
            return False

        detector_prefix = detector_prefixes[0]

        more_wanted_header_items = [
            "DETECTOR_DIMENSIONS",
            "DETECTOR_SIZE",
            "DETECTOR_VECTORS",
            "GONIO_NAMES",
            "GONIO_UNITS",
            "GONIO_VALUES",
            "GONIO_VECTORS",
        ]

        if any(
            f"{detector_prefix}{item}" not in header
            for item in more_wanted_header_items
        ):
            return False

        det_desc = "%sDETECTOR_DESCRIPTION" % detector_prefix
        if "CMOS-1" in header.get(det_desc, ""):
            return True

        return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment."""

        super().__init__(image_file, **kwargs)
        detector_prefixes = self._header_dictionary["DETECTOR_NAMES"].split()
        self._prefix = detector_prefixes[0]

    def _start(self):
        super()._start()
        self._header_size = int(self._header_dictionary["HEADER_BYTES"])

    def _goniometer(self):
        axis = tuple(map(float, self._header_dictionary["ROTATION_VECTOR"].split()))

        return self._goniometer_factory.known_axis(axis)

    def _detector(self):
        detector_name = self._header_dictionary["DETECTOR_NAMES"].split()[0].strip()

        detector_axes = self.get_detector_axes(detector_name)
        detector_fast = matrix.col(tuple(detector_axes[:3]))
        detector_slow = matrix.col(tuple(detector_axes[3:]))

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

        overload = int(float(self._header_dictionary["SATURATED_VALUE"]))
        underload = 0

        return self._detector_factory.complex(
            "CCD",
            detector_origin.elems,
            detector_fast.elems,
            detector_slow.elems,
            pixel_size,
            image_size,
            (underload, overload),
        )

    def _beam(self):
        beam_direction = self.get_beam_direction()
        p_fraction, p_plane = self.get_beam_polarization()

        wavelength = float(self._header_dictionary["WAVELENGTH"])

        return self._beam_factory.complex(
            beam_direction, p_fraction, p_plane, wavelength
        )

    def _scan(self):
        rotation = self.get_rotation()

        time_record = self._header_dictionary["DATE"]
        date_record, ms = time_record.split(".")
        epoch = calendar.timegm(time.strptime(date_record, "%a %b %d %Y %H:%M:%S"))
        epoch += 0.001 * int(ms)

        exposure_time = rotation[3]
        osc_start = rotation[0]
        osc_range = rotation[2]

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        # currently have no non-little-endian machines...
        assert len(self.get_detector()) == 1
        image_size = self.get_detector()[0].get_image_size()
        with self.open_file(self._image_file) as fh:
            fh.seek(self._header_size)
            raw_data = read_uint16(streambuf(fh), int(image_size[0] * image_size[1]))
        raw_data.reshape(flex.grid(image_size[1], image_size[0]))

        return raw_data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVCMOS1.understand(arg))
