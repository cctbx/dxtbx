from __future__ import absolute_import, division, print_function

import sys

import h5py

from scitbx.array_family import flex

from dxtbx.format.FormatHDF5 import FormatHDF5


class FormatRawData(FormatHDF5):
    @staticmethod
    def understand(image_file):
        with h5py.File(image_file, "r") as h5_handle:
            return len(h5_handle) == 1 and "data" in h5_handle

    def _start(self):
        self._h5_handle = h5py.File(self.get_image_file(), "r")

    def get_raw_data(self):
        data = self._h5_handle["data"]
        return flex.int(data[:, :])

    def get_num_images(self):
        return 1


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatRawData.understand(arg))
