"""
An implementation of the SMV image reader for Rigaku Eiger.
Be aware: this is completely unrelated to the HDF5 Eiger format.
"""

import sys
import time

from boost_adaptbx.boost.python import streambuf
from cctbx.eltbx import attenuation_coefficient
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.ext import read_int32
from dxtbx.format.FormatSMVRigaku import FormatSMVRigaku
from dxtbx.model import ParallaxCorrectedPxMmStrategy


class FormatSMVRigakuEiger(FormatSMVRigaku):
    """A class for reading SMV format Rigaku Eiger images."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like a Rigaku d*TREK SMVRigaku format image,
        i.e. we can make sense of it. Essentially that will be if it contains
        all of the keys we are looking for."""

        size, header = FormatSMVRigaku.get_smv_header(image_file)

        if "DETECTOR_NAMES" not in header:
            return False

        name = header["DETECTOR_NAMES"]

        if "%sDETECTOR_DESCRIPTION" % name not in header:
            return False

        return "Eiger" in header["%sDETECTOR_DESCRIPTION" % name]

    def _goniometer(self):
        """Initialize the structure for the goniometer."""

        values = [
            float(e) for e in self._header_dictionary["CRYSTAL_GONIO_VALUES"].split()
        ]
        names = [
            e.strip() for e in self._header_dictionary["CRYSTAL_GONIO_NAMES"].split()
        ]
        units = [
            e.strip() for e in self._header_dictionary["CRYSTAL_GONIO_UNITS"].split()
        ]
        axis_elts = [
            float(e) for e in self._header_dictionary["CRYSTAL_GONIO_VECTORS"].split()
        ]
        axes = [matrix.col(axis_elts[3 * j : 3 * (j + 1)]) for j in range(len(units))]
        scan_axis = self._header_dictionary["ROTATION_AXIS_NAME"].strip()

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

        gonio = self._goniometer_factory.make_multi_axis_goniometer(
            axes, values, names, scan_axis
        )

        # The calculated rotation axis is also recorded in the header. We could
        # use this to check that the goniometer is as expected
        rot_axis = tuple(map(float, self._header_dictionary["ROTATION_VECTOR"].split()))
        for e1, e2 in zip(rot_axis, gonio.get_rotation_axis()):
            assert abs(e1 - e2) < 1e-6

        return gonio

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
        underload = -1

        # Unfortunately thickness and material are not stored in the header. Set
        # to sensible defaults.
        t0 = 0.450
        material = "Si"

        table = attenuation_coefficient.get_table(material)
        wavelength = float(self._header_dictionary["SCAN_WAVELENGTH"])
        mu = table.mu_at_angstrom(wavelength) / 10.0

        detector = self._detector_factory.complex(
            "PAD",
            detector_origin.elems,
            detector_fast.elems,
            detector_slow.elems,
            pixel_size,
            image_size,
            (underload, overload),
            px_mm=ParallaxCorrectedPxMmStrategy(mu, t0),
            mu=mu,
        )

        for panel in detector:
            panel.set_thickness(t0)
            panel.set_material(material)

        return detector

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
        print(FormatSMVRigakuEiger.understand(arg))
