from __future__ import annotations

import ast
import optparse
import os
import sys
from urllib.request import urlretrieve

import procrunner
import py

import dxtbx.util


def find_format_classes(directory, base_python_path="dxtbx.format"):
    format_classes = []
    for name in directory.listdir("Format*.py"):
        content = name.read()
        try:
            parsetree = ast.parse(content)
        except SyntaxError:
            print("  *** Could not parse %s" % name.strpath)
            continue
        for top_level_def in parsetree.body:
            if not isinstance(top_level_def, ast.ClassDef):
                continue
            base_names = [
                baseclass.id
                for baseclass in top_level_def.bases
                if isinstance(baseclass, ast.Name)
            ]
            if any(n.startswith("Format") for n in base_names):
                classname = top_level_def.name
                format_classes.append(
                    "{classname}:{baseclasses} = {base_python_path}.{modulename}:{classname}".format(
                        baseclasses=",".join(base_names),
                        base_python_path=base_python_path,
                        classname=classname,
                        modulename=name.purebasename,
                    )
                )
                print("  found", classname, "based on", str(base_names))
    return format_classes


def install_package(home_location, format_classes):
    setup = (
        """
from __future__ import absolute_import, division, print_function

import os
import sys
from setuptools import setup

setup(
    author="Markus Gerstel",
    author_email="scientificsoftware@diamond.ac.uk",
    entry_points={
        "dxtbx.format": """
        + repr(format_classes)
        + """,
    },
    include_package_data=True,
    license="BSD license",
    name="dxtbx_custom",
    packages=["dxtbx_custom"],
    version="1.0.0",
    zip_safe=False,
)
"""
    )

    (home_location / "setup.py").write(setup)
    (home_location / "__init__.py").ensure()
    if not (home_location / "dxtbx_custom").check():
        (home_location / "dxtbx_custom").mksymlinkto(home_location)
    procrunner.run(
        ["libtbx.python", home_location / "setup.py", "develop", "--user"],
        working_directory=home_location,
        print_stdout=False,
    ).check_returncode()


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
    parser = optparse.OptionParser(
        usage="dxtbx.install_format (--user | --global) [/path/to/format/class.py] [URL]",
        description=(
            "Updates the dxtbx format class registry and installs format classes "
            "by copying/downloading them into the format class directory. "
            "The command must be re-run whenever format classes are added or removed "
            "manually."
        ),
    )
    parser.add_option("-?", action="help", help=optparse.SUPPRESS_HELP)
    parser.add_option(
        "-g",
        "--global",
        action="store_true",
        dest="glob",
        default=False,
        help=(
            "Install format classes globally (requires write access to base python, "
            "affects everyone using the installation, files go to /build/dxtbx/formats)"
        ),
    )
    parser.add_option(
        "-u",
        "--user",
        action="store_true",
        dest="user",
        default=False,
        help=(
            "Install format classes for current user only (requires write access to home directory, "
            "affects all python installations for the current user, files go to ~/.dxtbx)"
        ),
    )
    options, args = parser.parse_args(args)

    if options.glob:
        import libtbx.load_env

        home_location = py.path.local(abs(libtbx.env.build_path)) / "dxtbx" / "formats"
    elif options.user:
        home_location = py.path.local(os.path.expanduser("~")) / ".dxtbx"
    else:
        parser.print_help()
        print("\nYou must specify --global or --user")
        sys.exit(1 if args else 0)
    home_location.ensure(dir=True)

    for fc in args:
        try:
            local_file = py.path.local(fc)
        except Exception:
            local_file = False
        if local_file and local_file.check(file=1):
            home_location_copy = home_location / local_file.basename
            if local_file == home_location_copy:
                continue
            local_file.copy(home_location_copy)
            print("Copied", local_file.strpath, "to", home_location_copy.strpath)
            continue
        # Download the file from `url` and save it locally under `file_name`:
        if "/" in fc:
            local_file = home_location / fc.split("/")[-1]
            if local_file.ext != ".py":
                local_file = local_file.new(basename=local_file.basename + ".py")
            try:
                urlretrieve(fc, local_file.strpath)
            except Exception:
                local_file = False
            if local_file and local_file.check(file=1):
                print("Downloaded", fc, "to", local_file.strpath)
                continue
        sys.exit("Could not understand " + repr(fc))

    format_classes = find_format_classes(home_location, base_python_path="dxtbx_custom")
    install_package(home_location, format_classes)
    print("\nFormat classes installed, format class registry updated.")


if __name__ == "__main__":
    run()
