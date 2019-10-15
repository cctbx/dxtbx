"""
Print out the contents of the dxtbx understanding of a bunch of images to
an example XDS.INP file. This should illustrate the usage of the dxtbx
classes.
"""

from __future__ import absolute_import, division, print_function

import sys

from dxtbx.imageset import ImageSetFactory
from dxtbx.serialize import load, xds


def run(file_names):
    if len(file_names) == 1 and file_names[0].endswith("json"):
        try:
            datablock = load.datablock(file_names[0])
            assert len(datablock) == 1
            sequence = datablock[0].extract_sequences()[0]
        except ValueError as e:
            if str(e) == '"__id__" does not equal "imageset"':
                experiments = load.experiment_list(file_names[0])
                assert len(experiments) == 1
                sequence = experiments[0].imageset
            else:
                raise
    else:
        sequence = ImageSetFactory.new(file_names)[0]
    xsx = xds.to_xds(sequence)
    print(xsx.XDS_INP())


if __name__ == "__main__":
    run(sys.argv[1:])
