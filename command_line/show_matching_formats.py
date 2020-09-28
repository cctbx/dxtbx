from __future__ import absolute_import, division, print_function

import os
import sys

import dxtbx.format.Registry

dag = dxtbx.format.Registry.get_format_class_dag()


def recurse(parentformat, filename):
    for subformat in dag.get(parentformat, []):
        understood = dxtbx.format.Registry.get_format_class_for(subformat).understand(
            filename
        )
        print("%s: %s" % (subformat, understood))
        if understood:
            recurse(subformat, filename)


def show_matching_formats(files):
    for filename in files:
        print("\n=== %s ===" % filename)
        if os.path.exists(filename):
            recurse("Format", filename)
        else:
            print("File not found.")


def run(args=None):
    show_matching_formats(args or sys.argv[1:])


if __name__ == "__main__":
    run()
