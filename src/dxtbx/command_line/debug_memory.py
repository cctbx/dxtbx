# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.debug_memory

from __future__ import annotations

import argparse
import resource

import dxtbx
import dxtbx.util


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
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
