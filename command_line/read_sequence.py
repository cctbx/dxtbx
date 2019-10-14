# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.read_sequence

"""Tool to benchmark overall time cost for simply reading data"""

from __future__ import absolute_import, division, print_function

import sys
import time

from dxtbx.imageset import ImageSetFactory


def read_sequence(list_of_images):

    sequences = ImageSetFactory.new(list_of_images)

    for sequence in sequences:
        print(sequence.get_detector())
        print(sequence.get_scan())

        indices = sequence.indices()

        t0 = time.time()
        for i in indices:
            sequence.get_raw_data(i)
        t1 = time.time()

        print("Reading %d frames took %.2fs" % (len(indices), t1 - t0))


if __name__ == "__main__":
    read_sequence(sys.argv[1:])
