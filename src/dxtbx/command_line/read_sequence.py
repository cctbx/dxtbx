# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.read_sequence

"""Tool to benchmark overall time cost for simply reading data"""

from __future__ import annotations

import argparse
import time

import dxtbx.util
from dxtbx.imageset import ImageSetFactory


def read_sequence(images: list[str]):
    sequences = ImageSetFactory.new(images)

    for sequence in sequences:
        print(sequence.get_detector())
        print(sequence.get_scan())

        indices = sequence.indices()

        t0 = time.time()
        for i in indices:
            sequence.get_raw_data(i)
        t1 = time.time()

        print(f"Reading {len(indices)} frames took {t1 - t0:.2f}s")


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
    parser = argparse.ArgumentParser(
        description="Benchmark the time to read a set of images"
    )
    parser.add_argument("images", metavar="IMAGE", help="Images to read", nargs="+")
    options = parser.parse_args(args)
    return read_sequence(options.images)


if __name__ == "__main__":
    run()
