from __future__ import annotations

import ast
import optparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.request import urlretrieve

import dxtbx.util


def find_format_classes(
    directory: Path, base_python_path: str = "dxtbx.format"
) -> list[str]:
    format_classes = []
    for name in directory.glob("Format*.py"):
        content = name.read_bytes()
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
                print(f"  found {classname} based on {base_names}")
    return format_classes


def _ensure_cert_file():
    """
    Validate that the SSL_CERT_FILE environment variable is correct.

    On MacOS .pkg installations, the SSL_CERT_FILE path is incorrect.
    This means that the following GET fails. Check for this, and work
    out where certifi actually is.
    """
    cert_file = os.getenv("SSL_CERT_FILE")
    if cert_file and not Path(cert_file).resolve().is_file():
        try:
            import certifi

            os.environ["SSL_CERT_FILE"] = certifi.where()
        except ModuleNotFoundError:
            # Better to try default instead of a nonexistent file
            del os.environ["SSL_CERT_FILE"]


def install_package(home_location: Path, format_classes: list[str]):
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

    (home_location / "setup.py").write_text(setup)
    (home_location / "__init__.py").touch()
    custom_folder = home_location / "dxtbx_custom"
    if not custom_folder.exists():
        os.symlink(home_location, custom_folder)
    subprocess.run(
        ["libtbx.python", str(home_location / "setup.py"), "develop", "--user"],
        cwd=home_location,
        capture_output=True,
        check=True,
    )


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

        home_location = Path(abs(libtbx.env.build_path)) / "dxtbx" / "formats"
    elif options.user:
        home_location = Path.home() / ".dxtbx"
    else:
        parser.print_help()
        print("\nYou must specify --global or --user")
        sys.exit(1 if args else 0)
    home_location.mkdir(exist_ok=True, parents=True)

    for fc in args:
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
            _ensure_cert_file()
            local_file = home_location / (Path(fc.split("/")[-1]).stem + ".py")
            try:
                urlretrieve(fc, str(local_file))
            except Exception as e:
                print("Could not load from URL:", e)
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
