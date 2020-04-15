from __future__ import absolute_import, division, print_function

import sys

import dxtbx.format.Registry

dag = dxtbx.format.Registry.get_format_class_dag()


def print_class(class_name, filename=None, depth=1):
    if filename is None:
        print("% 5d" % depth, "  " * depth, class_name)
    else:
        format_class = dxtbx.format.Registry.get_format_class_for(class_name)
        if format_class.understand(filename):
            print("% 5d" % depth, "  " * depth, class_name)
        else:
            return
    for child in dag.get(class_name, []):
        print_class(child, filename, depth + 1)


def show_registry(filename=None):
    if filename is None:
        extrabit = ""
    else:
        extrabit = " that understand image %s" % filename
    print(
        "Showing hierarchy of classes in the dxtbx registry%s. The root classes are shown with depth of 1, and subclasses are shown indented and with a higher depth number."
        % extrabit
    )
    print()
    print("Depth  Class name")
    print("    0  Format")

    for class_name in dag["FormatFile"]:
        print_class(class_name, filename)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        show_registry(sys.argv[1])
    else:
        show_registry()
