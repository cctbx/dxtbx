# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.debug_memory

import argparse
import resource
from builtins import range

import dxtbx


def run(args=None):
    parser = argparse.ArgumentParser(
        description="Test memory usage by repeatedly loading an image"
    )
    parser.add_argument("image_frame", help="An image to reload repeatedly")
    options = parser.parse_args(args)

    powers = [2 ** n for n in range(20)]

    for j in range(powers[-1] + 1):
        dxtbx.load(options.image_frame)
        mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if j in powers:
            print(j, mem)


if __name__ == "__main__":
    run()
