from __future__ import absolute_import, division, print_function

import h5py

from dxtbx.format.FormatNexus import FormatNexus


class FormatNexusTimepix2M(FormatNexus):
    @staticmethod
    def understand(image_file):
        # Get the file handle
        with h5py.File(image_file, "r") as handle:
            if "/entry/instrument/detector/module/data_size" in handle:
                size = handle["/entry/instrument/detector/module/data_size"]
                if tuple(size[()]) == (1147, 2069):
                    return True

        return False

    def __init__(self, image_file, **kwargs):

        super(FormatNexusTimepix2M, self).__init__(image_file, **kwargs)

    def _detector(self):
        """return a detector with additional masking"""

        detector = super(FormatNexusTimepix2M, self)._detector()

        # add masks for intermediate regions:

        module_spacing = 117
        chip_spacing = 3

        mask_width = 5

        nfast = 8 * 256 + 7 * chip_spacing
        nslow = 2 * 256 + chip_spacing

        for module in 0, 1:
            # horizontal masks
            for j in range(1):
                nn = module * (nslow + module_spacing)
                detector[0].add_mask(0, nn + 254, nfast, nn + 254 + mask_width)
            # vertical masks
            for j in range(1, 8):
                nn = module * (nslow + module_spacing)
                mm = j * 259
                detector[0].add_mask(mm - 5, nn, mm, nn + nslow)

        # between modules
        detector[0].add_mask(0, nslow, nfast, nslow + module_spacing)

        return detector
