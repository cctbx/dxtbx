# LIBTBX_SET_DISPATCHER_NAME dxtbx.show_mask_info

import sys

from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.util import show_mask_info


def run(args):
    filenames = args or sys.argv[1:]
    try:
        el = ExperimentListFactory.from_filenames(filenames)
    except FileNotFoundError as e:
        sys.exit(str(e))

    show_mask_info(el)


if __name__ == "__main__":
    run()
