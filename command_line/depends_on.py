import argparse

import h5py

import dxtbx.util

sample = None


def depends_on(in_name):
    f = h5py.File(in_name, "r")

    depends = {}

    global sample

    def finder(thing, path):
        if hasattr(thing, "attrs"):
            if "depends_on" in thing.attrs:
                depends[path] = thing.attrs["depends_on"]
            if thing.attrs.get("NX_class", None) == "NXsample":
                global sample
                sample = path
        if hasattr(thing, "keys"):
            for k in thing:
                try:
                    finder(thing[k], path=f"{path}/{k}")
                except (IOError, TypeError, ValueError, KeyError):
                    pass

    # clean up hierarchy to just have sample stuff

    finder(f, path="")
    delete = []

    for d in sorted(depends):
        if "/entry/sample/transformations" not in d:
            delete.append(d)

    for d in delete:
        del depends[d]

    # invert for printing

    inverted = {}

    for d in depends:
        t = depends[d]
        if t in inverted:
            print(t, inverted[t], d)
        inverted[t] = d

    print("Dependency hierarchy in file:")
    at = "."
    depth = 0
    while at in inverted:
        print(f"{' ' * depth}+ {at}")
        at = inverted[at]
        depth += 2
    print(f"{' ' * depth}+ {at}")
    print("")

    print(f"Sample at {sample} depends on:")
    if "depends_on" in f[sample]:
        print(f[sample]["depends_on"][()])
    elif hasattr(f[sample], "attrs") and "depends_on" in f[sample].attrs:
        print(f[sample].attrs["depends_on"])
    else:
        print(f"{sample} -> depends_on not found")

    f.close()


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
    parser = argparse.ArgumentParser(
        description="Print depends_on hierarchy for Nexus files"
    )
    parser.add_argument("filename", help="The nexus file")
    opts = parser.parse_args(args)
    depends_on(opts.filename)


if __name__ == "__main__":
    run()
