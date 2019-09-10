from __future__ import absolute_import, division, print_function

from iotbx.detectors.dtrek import DTREKImage

from dxtbx.format.FormatSMVRigakuA200 import FormatSMVRigakuA200


class FormatSMVRigakuA200SPring8BL26B1(FormatSMVRigakuA200):
    """A class for reading SMV format Rigaku A200 images written by a detector
    usually at SPring-8 BL26B1, which requires a reversed axis goniometer"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like a Rigaku A200 SMV format image,
        i.e. we can make sense of it. Essentially that will be if it contains
        all of the keys we are looking for."""

        size, header = FormatSMVRigakuA200SPring8BL26B1.get_smv_header(image_file)

        detector_prefix = header["DETECTOR_NAMES"].split()[0].strip()

        test = "%s%s" % (detector_prefix, "DETECTOR_IDENTIFICATION")
        if header.get(test) != "MSC_REIT_A200_SN09250183":
            return False

        test = "%s%s" % (detector_prefix, "DETECTOR_DESCRIPTION")
        if not header.get(test).startswith("A200"):
            return False

        test = "%s%s" % (detector_prefix, "GONIO_DESCRIPTION")
        if not header.get(test).startswith("SPring-8"):
            return False

        return True

    def detectorbase_start(self):
        self.detectorbase = DTREKImage(self._image_file)
        self.detectorbase.readHeader()

    def _goniometer(self):
        """Initialize the structure for the goniometer. Although ROTATION_VECTOR
        is written into the header, this might (?) always be 1. 1. 0. for this
        detector, and we want a reversed axis goniometer"""

        # axis = tuple(map(float, self._header_dictionary[
        #    'ROTATION_VECTOR'].split()))

        return self._goniometer_factory.single_axis_reverse()
