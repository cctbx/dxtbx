# LIBTBX_SET_DISPATCHER_NAME dxtbx.show_mask_info

import argparse
import sys

from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.util import show_mask_info


def run(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", metavar="IMAGE", nargs="+")
    options = parser.parse_args(args)
    try:
        el = ExperimentListFactory.from_filenames(options.filenames)
    except FileNotFoundError as e:
        sys.exit(str(e))

    show_mask_info(el)


if __name__ == "__main__":
    run()
