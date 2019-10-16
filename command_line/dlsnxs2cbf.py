from __future__ import absolute_import, division, print_function

import sys

import dxtbx.util.dlsnxs2cbf

if __name__ == "__main__":
    dxtbx.util.dlsnxs2cbf.make_cbf(sys.argv[1], sys.argv[2])
