from __future__ import annotations

import argparse
import ast
import shutil
import sys
from pathlib import Path
from urllib.request import urlretrieve

import procrunner

import dxtbx.util


def find_format_classes(directory: Path, base_python_path="dxtbx.format"):
    format_classes = []
    for name in directory.glob("Format*.py"):
        content = name.read_text()
        try:
            parsetree = ast.parse(content)
        except SyntaxError:
            print(f"  *** Could not parse {name}")
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
                        modulename=name.stem,
                    )
                )
                print("  found", classname, "based on", str(base_names))
    return format_classes


def install_package(home_location: Path, format_classes):
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

    home_location.joinpath("setup.py").write_text(setup)
    home_location.joinpath("__init__.py").touch()
    if not (dxtbx_custom := home_location / "dxtbx_custom").exists():
        dxtbx_custom.symlink_to(home_location, target_is_directory=True)
    procrunner.run(
        ["libtbx.python", home_location / "setup.py", "develop", "--user"],
        working_directory=home_location,
        print_stdout=False,
    ).check_returncode()


def run(args=None):
    # usage="dxtbx.install_format (--user | --global) [/path/to/format/class.py] [URL]",
    dxtbx.util.encode_output_as_utf8()
    parser = argparse.ArgumentParser(
        description=(
            "Updates the dxtbx format class registry and installs format classes "
            "by copying/downloading them into the format class directory. "
            "The command must be re-run whenever format classes are added or removed "
            "manually."
        ),
    )
    parser.add_argument("-?", action="help", help=argparse.SUPPRESS)
    parser.add_argument(
        "sources",
        metavar="SOURCE",
        nargs="+",
        help="The path to the format class file (/path/to/format/class.py) to import, or URL to download from",
    )
    destination_group = parser.add_mutually_exclusive_group(required=True)
    destination_group.add_argument(
        "-g",
        "--global",
        action="store_true",
        dest="write_global",
        help=(
            "Install format classes globally (requires write access to base python, "
            "affects everyone using the installation, files go to /build/dxtbx/formats)"
        ),
    )
    destination_group.add_argument(
        "-u",
        "--user",
        action="store_true",
        dest="write_user",
        help=(
            "Install format classes for current user only (requires write access to home directory, "
            "affects all python installations for the current user, files go to ~/.dxtbx)"
        ),
    )
    args = parser.parse_args(args)

    if args.write_global:
        import libtbx.load_env

        home_location = Path(abs(libtbx.env.build_path)) / "dxtbx" / "formats"
    elif args.write_user:
        home_location = Path.home() / ".dxtbx"

    home_location.mkdir(exist_ok=True)

    for fc in args.sources:
        try:
            local_file = Path(fc)
        except Exception:
            local_file = False
        if local_file and local_file.is_file():
            home_location_copy = home_location / local_file.name
            if local_file == home_location_copy:
                continue
            shutil.copy(local_file, home_location_copy)
            print(f"Copied {local_file} to {home_location_copy}")
            continue
        # Download the file from `url` and save it locally under `file_name`:
        if "/" in fc:
            local_file = home_location / fc.split("/")[-1]
            if local_file.suffix != ".py":
                local_file = local_file.parent / (local_file.name + ".py")
            try:
                urlretrieve(fc, str(local_file))
            except Exception:
                local_file = False
            if local_file and local_file.is_file():
                print(f"Downloaded {fc} to {local_file}")
                continue
        sys.exit("Could not understand " + repr(fc))

    format_classes = find_format_classes(home_location, base_python_path="dxtbx_custom")
    install_package(home_location, format_classes)
    print("\nFormat classes installed, format class registry updated.")


if __name__ == "__main__":
    run()
