from __future__ import annotations

import argparse

import dxtbx.util
from dxtbx import load


def saturation(image_file):
    i = load(image_file)
    d = i.get_detector()
    raw_data = i.get_raw_data()
    if not isinstance(raw_data, tuple):
        raw_data = (raw_data,)
    if i.get_scan() is None:
        return (
            0,
            max(
                max(raw_data[pid]) / detector.get_trusted_range()[1]
                for pid, detector in enumerate(d)
            ),
        )
    else:
        return (
            i.get_scan().get_image_range()[0],
            max(
                max(raw_data[pid]) / detector.get_trusted_range()[1]
                for pid, detector in enumerate(d)
            ),
        )


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
    parser = argparse.ArgumentParser()
    parser.add_argument("images", metavar="IMAGE", help="Image files", nargs="+")
    options = parser.parse_args(args)

    for image_file in options.images:
        i, s = saturation(image_file)
        print(f"{i:6d} {s:.6f}")


if __name__ == "__main__":
    run()
