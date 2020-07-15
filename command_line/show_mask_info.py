# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.show_mask_info

import sys

from dxtbx.model.experiment_list import ExperimentListFactory


def show_mask_info():
    try:
        el = ExperimentListFactory.from_filenames(sys.argv[1:])
    except FileNotFoundError as e:
        sys.exit(str(e))

    for i in el.imagesets():
        d = i.get_detector()
        m = i.get_mask(0)
        print(f"---- ----")
        print(d)
        for j, _m in enumerate(m):
            print(f"Module {j} has {_m.count(False)} masked pixels of {_m.size()}")


if __name__ == "__main__":
    show_mask_info()
