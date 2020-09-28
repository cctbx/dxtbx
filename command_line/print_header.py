from __future__ import absolute_import, division, print_function

import sys

from scitbx.array_family import flex

import dxtbx.format.Registry
from dxtbx.format.FormatMultiImage import FormatMultiImage


def print_header(args=None):
    args = args or sys.argv[1:]
    # this will do the lookup for every frame - this is strictly not needed
    # if all frames are from the same instrument
    for arg in args:
        print("=== %s ===" % arg)
        format_class = dxtbx.format.Registry.get_format_class_for_file(arg)
        if not format_class:
            print("No format class found that can understand %s" % arg)
            continue
        print("Using header reader: %s" % format_class.__name__)
        i = format_class(arg)
        beam = i.get_beam()
        goniometer = i.get_goniometer()
        detector = i.get_detector()
        scan = i.get_scan()
        if beam is None:
            print("No beam model found")
        else:
            print(beam)
        if detector is None:
            print("No detector model found")
        else:
            print(detector)
        if goniometer is None:
            print("No goniometer model found")
        else:
            print(goniometer)
        if scan is None:
            print("No scan model found")
        else:
            print(scan)

        if not issubclass(format_class, FormatMultiImage):
            try:
                raw_data = i.get_raw_data()
                if not isinstance(raw_data, tuple):
                    raw_data = (raw_data,)
                d = [p.as_1d() for p in raw_data]
                print("Total Counts: %d" % sum([flex.sum(p.select(p >= 0)) for p in d]))
            except AttributeError:
                print("Could not read image data")


def run(args=None):
    print_header(args)


if __name__ == "__main__":
    run()
