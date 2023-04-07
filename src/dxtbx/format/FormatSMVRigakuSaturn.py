"""
An implementation of the SMV image reader for Rigaku Saturn images.
Inherits from FormatSMVRigaku.
"""

from __future__ import annotations

import sys
import time

from iotbx.detectors.saturn import SaturnImage
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.format.FormatSMVRigaku import FormatSMVRigaku


class FormatSMVRigakuSaturn(FormatSMVRigaku):
    """A class for reading SMV format Rigaku Saturn images, and correctly
    constructing a model for the experiment from this."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like a Rigaku Saturn SMV format image,
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
            "DTREK_DATE_TIME",
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

        descriptive_items = ["DETECTOR_IDENTIFICATION", "DETECTOR_DESCRIPTION"]

        for header_item in descriptive_items:
            test = f"{detector_prefix}{header_item}"
            if test in header and "saturn" in header[test].lower():
                return True

        return False

    def detectorbase_start(self):
        self.detectorbase = SaturnImage(self._image_file)
        self.detectorbase.open_file = self.open_file
        self.detectorbase.readHeader()

    def _goniometer(self):
        """Initialize the structure for the goniometer - this will need to
        correctly compose the axes given in the image header. In this case
        this is made rather straightforward as the image header has the
        calculated rotation axis stored in it. We work from the
        rest of the header and construct a goniometer model."""

        head_dict = self._header_dictionary

        names = [e.strip() for e in head_dict["CRYSTAL_GONIO_NAMES"].split()]
        values = [float(e) for e in head_dict["CRYSTAL_GONIO_VALUES"].split()]
        units = [e.strip() for e in head_dict["CRYSTAL_GONIO_UNITS"].split()]
        axis_elts = [float(e) for e in head_dict["CRYSTAL_GONIO_VECTORS"].split()]

        rot_axis = tuple(map(float, head_dict["ROTATION_VECTOR"].split()))
        scan_axis = head_dict["ROTATION_AXIS_NAME"].strip()
        axes = [matrix.col(axis_elts[3 * j : 3 * (j + 1)]) for j in range(len(units))]

        # Take only elements that have corresponding units of 'deg' (which is
        # probably all of them).
        filt = [e == "deg" for e in units]
        values = [e for e, f in zip(values, filt) if f]
        names = [e for e, f in zip(names, filt) if f]
        axes = [e for e, f in zip(axes, filt) if f]

        # Multi-axis gonio requires axes in order as viewed from crystal to gonio
        # base. Assume the SMV header records them in reverse order.

        axes = flex.vec3_double(reversed(axes))
        names = flex.std_string(reversed(names))
        values = flex.double(reversed(values))
        scan_axis = flex.first_index(names, scan_axis)
        assert scan_axis is not None
        gonio = self._goniometer_factory.make_multi_axis_goniometer(
            axes, values, names, scan_axis
        )

        # The calculated rotation axis is also recorded in the header. We
        # use this to check that the goniometer is as expected
        for e1, e2 in zip(rot_axis, gonio.get_rotation_axis()):
            assert abs(e1 - e2) < 1e-6

        return gonio

    def _detector(self):
        """Return a model for the detector, allowing for two-theta offsets
        and the detector position. This will be rather more complex..."""

        detector_name = self._header_dictionary["DETECTOR_NAMES"].split()[0].strip()

        detector_axes = self.get_detector_axes(detector_name)
        detector_fast = matrix.col(tuple(detector_axes[:3]))
        detector_slow = matrix.col(tuple(detector_axes[3:]))

        # apply spatial distortion
        distortion = self.get_distortion(detector_name)
        detector_fast, detector_slow = (
            distortion[0] * detector_fast + distortion[1] * detector_slow,
            distortion[2] * detector_fast + distortion[3] * detector_slow,
        )

        beam_pixels = self.get_beam_pixels(detector_name)
        pixel_size = self.get_pixel_size(detector_name)
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

        image_size = self.get_image_size(detector_name)

        rotations = []
        translations = []

        for j, unit in enumerate(gonio_units):
            axis = matrix.col(gonio_axes[3 * j : 3 * (j + 1)])
            value = gonio_values[j]

            if unit == "deg":
                rotations.append(
                    axis.axis_and_angle_as_r3_rotation_matrix(value, deg=True)
                )
                translations.append(matrix.col((0.0, 0.0, 0.0)))
            elif unit == "mm":
                rotations.append(
                    matrix.sqr((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
                )
                translations.append(value * axis)
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
        epoch = time.strptime(
            self._header_dictionary["DTREK_DATE_TIME"], "%d-%b-%Y %H:%M:%S"
        )
        return self._create_single_SVM_scan(epoch, local_time=False)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVRigakuSaturn.understand(arg))
