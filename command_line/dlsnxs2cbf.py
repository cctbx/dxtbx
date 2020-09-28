from __future__ import absolute_import, division, print_function

import sys

import dxtbx.util.dlsnxs2cbf


def run(args=None):
    args = args or sys.argv[1:]
    dxtbx.util.dlsnxs2cbf.make_cbf(args[0], args[1])


if __name__ == "__main__":
    run()
