# LIBTBX_SET_DISPATCHER_NAME dev.dxtbx.show_mask_info

import sys

from dxtbx.model.experiment_list import ExperimentListFactory


def show_mask_info():
    el = ExperimentListFactory.from_filenames(sys.argv[1:])

    for j, e in enumerate(el):
        i = e.imageset
        d = e.detector
        m = i.get_mask(0)
        print(f"---- Experiment {j} ----")
        print(d)
        for i, _m in enumerate(m):
            print(f"Module {i} has {_m.count(False)} masked pixels of {_m.size()}")


if __name__ == "__main__":
    show_mask_info()
