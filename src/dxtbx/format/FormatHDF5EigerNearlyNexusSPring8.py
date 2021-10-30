import sys

import h5py

from dxtbx.format.FormatHDF5EigerNearlyNexus import FormatHDF5EigerNearlyNexus


class FormatHDF5EigerNearlyNexusSPring8(FormatHDF5EigerNearlyNexus):
    @staticmethod
    def understand(image_file):
        is_nexus = FormatHDF5EigerNearlyNexus.understand(image_file)
        if not is_nexus:
            return False

        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            if "/entry/instrument/detector/detector_number" in handle and handle[
                "/entry/instrument/detector/detector_number"
            ][()] in [b"E-32-0114", b"E-32-0112", b"E-18-0103"]:
                return True
        return False

    def _start(self):
        super()._start()
        # invert the rotation axis
        rotation_axis = self._goniometer_model.get_rotation_axis_datum()
        self._goniometer_model.set_rotation_axis_datum([-x for x in rotation_axis])


if __name__ == "__main__":
    print(FormatHDF5EigerNearlyNexusSPring8.understand(sys.argv[1]))
    print(FormatHDF5EigerNearlyNexusSPring8(sys.argv[1]).get_goniometer())
