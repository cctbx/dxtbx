import argparse
import os

import dxtbx.format.Registry

dag = dxtbx.format.Registry.get_format_class_dag()


def recurse(parentformat, filename):
    for subformat in dag.get(parentformat, []):
        understood = dxtbx.format.Registry.get_format_class_for(subformat).understand(
            filename
        )
        print(f"{subformat}: {understood}")
        if understood:
            recurse(subformat, filename)


def show_matching_formats(files):
    for filename in files:
        print(f"\n=== {filename} ===")
        if os.path.exists(filename):
            recurse("Format", filename)
        else:
            print("File not found.")


def run(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", metavar="IMAGE", nargs="+")
    options = parser.parse_args(args)
    show_matching_formats(options.filenames)


if __name__ == "__main__":
    run()
