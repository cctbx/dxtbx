"""
An implementation of the SMV image reader for Rigaku images.
Inherits from FormatSMV.
"""

import calendar
import sys
import time

from dxtbx.format.FormatSMV import FormatSMV


class FormatSMVRigaku(FormatSMV):
    """A class for reading SMV format Rigaku images. 'Abstract' class."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like a Rigaku d*TREK SMV format image,
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

        if any(header_item not in header for header_item in wanted_header_items):
            return False

        # code around CMOS1 on ALS 4.2.2

        detector_prefixes = header["DETECTOR_NAMES"].split()

        if len(detector_prefixes) == 1:
            prefix = detector_prefixes[0]
            det_desc = "%sDETECTOR_DESCRIPTION" % prefix
            if "CMOS-1" in header.get(det_desc, ""):
                return False

        detector_prefixes = header["DETECTOR_NAMES"].split()
        try:
            detector_number = int(header["DETECTOR_NUMBER"].strip())
        except (KeyError, AttributeError, ValueError):
            return False

        if detector_number != len(detector_prefixes):
            return False

        return True

    def _goniometer(self):
        """Overload this method to read the image file however you like so
        long as the result is an goniometer."""

        raise NotImplementedError("overload me")

    def _detector(self):
        """Overload this method to read the image file however you like so
        long as the result is an detector."""

        raise NotImplementedError("overload me")

    def _beam(self):
        """Overload this method to read the image file however you like so
        long as the result is an beam."""

        raise NotImplementedError("overload me")

    def _scan(self):
        """Overload this method to read the image file however you like so
        long as the result is an scan."""

        raise NotImplementedError("overload me")

    def _create_single_SVM_scan(self, epoch_time_struct, local_time=True):
        """Return the scan information for this image."""

        rotation = self.get_rotation()

        if local_time:
            epoch = time.mktime(epoch_time_struct)
        else:
            epoch = calendar.timegm(epoch_time_struct)

        osc_start = rotation[0]
        osc_range = rotation[2]
        exposure_time = rotation[3]

        return self._scan_factory.single_file(
            self._image_file, exposure_time, osc_start, osc_range, epoch
        )

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
        flex array of integers.)"""

        assert len(self.get_detector()) == 1
        panel = self.get_detector()[0]
        image_size = panel.get_image_size()
        return self._get_endianic_raw_data(size=image_size)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatSMVRigaku.understand(arg))
