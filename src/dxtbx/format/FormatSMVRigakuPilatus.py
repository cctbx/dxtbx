"""
An implementation of the SMV image reader for Rigaku Pilatus images.
Inherits from FormatSMVRigaku.
"""

import sys
import time

from boost_adaptbx.boost.python import streambuf
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.ext import read_int32
from dxtbx.format.FormatSMVRigaku import FormatSMVRigaku


class FormatSMVRigakuPilatus(FormatSMVRigaku):
    """A class for reading SMV format Rigaku Pilatus images."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like a Rigaku d*TREK SMVRigaku format image,
        i.e. we can make sense of it. Essentially that will be if it contains
        all of the keys we are looking for."""

        size, header = FormatSMVRigaku.get_smv_header(image_file)

        if "DETECTOR_TYPE" not in header:
            return False

        if header["DETECTOR_TYPE"] not in ["Pilatus 200K", "Pilatus 300K"]:
            return False

        return True

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
        fast = matrix.col(tuple(detector_axes[:3]))
        slow = matrix.col(tuple(detector_axes[3:]))
        distortion = self.get_distortion(detector_name)

        # multiply through by the distortion to get the true detector fast, slow

        detector_fast, detector_slow = (
            distortion[0] * fast + distortion[1] * slow,
            distortion[2] * fast + distortion[3] * slow,
        )

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
            "PAD",
            detector_origin.elems,
            detector_fast.elems,
            detector_slow.elems,
            pixel_size,
            image_size,
            (underload, overload),
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
        epoch = time.strptime(
            self._header_dictionary["DTREK_DATE_TIME"], "%d-%b-%Y %H:%M:%S"
        )
        return self._create_single_SVM_scan(epoch)

    def get_raw_data(self):
        """Read the data - assuming it is streaming 4-byte unsigned ints following the
        header..."""

        assert len(self.get_detector()) == 1
        size = self.get_detector()[0].get_image_size()
        f = self.open_file(self._image_file)
        f.read(int(self._header_dictionary["HEADER_BYTES"]))
        raw_data = read_int32(streambuf(f), int(size[0] * size[1]))
        raw_data.reshape(flex.grid(size[1], size[0]))

        return raw_data


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVRigakuPilatus.understand(arg))
