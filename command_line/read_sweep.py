# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.read_sweep

"""Tool to benchmark overall time cost for simply reading data"""

from __future__ import absolute_import, division, print_function

import time


def read_sweep(list_of_images):

    from dxtbx.imageset import ImageSetFactory

    sweeps = ImageSetFactory.new(list_of_images)

    for sweep in sweeps:
        print(sweep.get_detector())
        print(sweep.get_scan())

        indices = sweep.indices()

        t0 = time.time()
        for i in indices:
            sweep.get_raw_data(i)
        t1 = time.time()

        print("Reading %d frames took %.2fs" % (len(indices), t1 - t0))


if __name__ == "__main__":
    import sys

    read_sweep(sys.argv[1:])
