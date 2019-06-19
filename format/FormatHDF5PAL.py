from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatHDF5 import FormatHDF5

"""
Format class for the PAL XFEL raw data format
"""

class FormatHDF5PAL(FormatHDF5):
    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatHDF5.__init__(self, image_file, **kwargs)

    @staticmethod
    def understand(image_file):
        import h5py

        h5_handle = h5py.File(image_file, "r")
        if len(h5_handle) != 1: return False

        for key in h5_handle:
            if not key.startswith('R'): return False
            try:
                run = int(key.lstrip('R'))
            except ValueError:
                return False
            for subkey in h5_handle[key]:
                if subkey not in ('header', 'scan_dat'): return False
                if subkey == 'scan_dat':
                    if 'raymx_header' not in h5_handle[key][subkey]: return False
        return True

    def _start(self):
        import h5py

        self._h5_handle = h5py.File(self.get_image_file(), "r")
        self._run = self._h5_handle.keys()[0]

        # currently hardcoded to Rayonix MX225HS
        self._detector_size = 225 # mm
        self._max_pixels = 5760
        frame_1 = self.get_raw_data(0)
        assert frame_1.focus()[0] == frame_1.focus()[1]
        self._binning = self._max_pixels // frame_1.focus()[0]

    def get_raw_data(self, index = None):
        from scitbx.array_family import flex
        import numpy as np
        if index is None: index = 0

        data = self._h5_handle[self._run]['scan_dat/raymx_data'][index]
        #return flex.int(int) # this crashes!
        #return flex.int(data.astype(np.int)) # this doesn't work! (data is read incorrectly)
        return flex.double(data.astype(np.float))

    def get_num_images(self):
        return len(self._h5_handle[self._run]['scan_dat/N'][()])

    def _detector(self):
        distance = self._h5_handle[self._run]['header/detector_0_distance'][()]
        image_size = self._max_pixels // self._binning
        pixel_size = self._detector_size / image_size
        beam_x = 0.5 * self._detector_size
        beam_y = 0.5 * self._detector_size
        trusted_range = (-1, 65534) # note one less than normal due to how bad pixels are encoded for this detector

        return self._detector_factory.simple(
            "CCD",
            distance,
            (beam_y, beam_x),
            "+x",
            "-y",
            (pixel_size, pixel_size),
            (image_size, image_size),
            trusted_range,
            [],
        )

    def _beam(self, index = None):
       if index is None: index = 0
       return self._beam_factory.simple(self._h5_handle[self._run]['scan_dat/photon_wavelength'][index])

if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        print(FormatHDF5PAL.understand(arg))
