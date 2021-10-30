# LIBTBX_SET_DISPATCHER_NAME dxtbx.show_mask_info

import argparse
import sys

import dxtbx.util
from dxtbx.model.experiment_list import ExperimentListFactory


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", metavar="IMAGE", nargs="+")
    options = parser.parse_args(args)
    try:
        el = ExperimentListFactory.from_filenames(options.filenames)
    except FileNotFoundError as e:
        sys.exit(str(e))

    dxtbx.util.show_mask_info(el)


if __name__ == "__main__":
    run()
