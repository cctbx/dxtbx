# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.print_beam_centre
"""
Print out the contents of the dxtbx understanding of a bunch of images to
an example XDS.INP file. This should illustrate the usage of the dxtbx
classes.
"""
from __future__ import absolute_import, division, print_function

import sys

from dxtbx.imageset import ImageSetFactory
from dxtbx.serialize import load


def run(file_names):
    if len(file_names) == 1 and file_names[0].endswith("json"):
        datablock = load.datablock(file_names[0])
        assert len(datablock) == 1
        sequence = datablock[0].extract_sequences()[0]
    else:
        sequence = ImageSetFactory.new(file_names)[0]
    detector = sequence.get_detector()
    beam = sequence.get_beam()
    print(detector.get_ray_intersection(beam.get_s0())[1])


if __name__ == "__main__":
    run(sys.argv[1:])
