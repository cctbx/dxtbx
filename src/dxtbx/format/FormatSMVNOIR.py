"""
An implementation of the SMV image reader for Rigaku Saturn images.
Inherits from FormatSMVRigaku.
"""


import sys
import time

from iotbx.detectors.noir import NoirImage
from scitbx import matrix

from dxtbx.format.FormatSMVRigaku import FormatSMVRigaku


class FormatSMVNOIR(FormatSMVRigaku):
    """A class for reading SMV format ALS 4.2.2 NOIR images, and correctly
    constructing a model for the experiment from this."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like a ALS 4.2.2 NOIR SMV format image,
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
            "NOIR1_CREATED",
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
        ]

        if any(
            f"{detector_prefix}{item}" not in header
            for item in more_wanted_header_items
        ):
            return False

        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
        proper model of the experiment. Easy from Rigaku Saturn images as
        they contain everything pretty much we need..."""

        super().__init__(image_file, **kwargs)
        self.detector_class = "NOIR1"
        self.detector = "adsc"

    def detectorbase_start(self):
        self.detectorbase = NoirImage(self._image_file)
        self.detectorbase.open_file = self.open_file
        self.detectorbase.readHeader()

    def _goniometer(self):
        """Initialize the structure for the goniometer - this will need to
        correctly compose the axes given in the image header. In this case
        this is made rather straightforward as the image header has the
        calculated rotation axis stored in it. We could work from the
        rest of the header and construct a goniometer model."""

        axis = tuple(map(float, self._header_dictionary["ROTATION_VECTOR"].split()))

        return self._goniometer_factory.known_axis(axis)

    def get_beam_pixels(self, detector_name):
        return [
            float(bp)
            for bp in self._header_dictionary[
                "%sSPATIAL_BEAM_POSITION" % detector_name
            ].split()[:2]
        ]

    def _detector(self):
        """Return a model for the detector, allowing for two-theta offsets
        and the detector position. This will be rather more complex..."""

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
        """Return a simple model for the beam."""

        beam_direction = self.get_beam_direction()
        p_fraction, p_plane = self.get_beam_polarization()

        wavelength = float(self._header_dictionary["SOURCE_WAVELENGTH"].split()[-1])

        return self._beam_factory.complex(
            beam_direction, p_fraction, p_plane, wavelength
        )

    def _scan(self):
        """Return the scan information for this image."""
        epoch = time.strptime(
            self._header_dictionary["NOIR1_CREATED"], "%m/%d/%y  %H:%M:%S"
        )
        return self._create_single_SVM_scan(epoch)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVNOIR.understand(arg))
