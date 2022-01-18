import argparse
from typing import Optional

import dxtbx.util
from dxtbx.format import Registry

dag = Registry.get_format_class_dag()


def print_class(class_name: str, filename: Optional[str] = None, depth=1):
    """Print a Format class name if it matches a file"""
    if filename is None or (
        Registry.get_format_class_for(class_name).understand(filename)
    ):
        print(f"{depth: 5} {'  ' * depth} {class_name}")
    else:
        return

    for child in dag.get(class_name, []):
        print_class(child, filename, depth + 1)


def show_registry(filename: Optional[str] = None):
    if filename:
        print(f"Format classes that understand {filename}:")
    else:
        print("Format classes in the dxtbx registry:")

    print("Depth  Class name")
    print("    0  Format")

    for class_name in dag["Format"]:
        print_class(class_name, filename)


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
    parser = argparse.ArgumentParser(
        description="Show hierarchy of dxtbx format classes"
    )
    parser.add_argument(
        "filename", help="Only show classes that understand this file", nargs="?"
    )
    opts = parser.parse_args(args)
    show_registry(opts.filename)


if __name__ == "__main__":
    run()
