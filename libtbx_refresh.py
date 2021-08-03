import os
import site
import subprocess
import sys
from pathlib import Path

import libtbx.load_env
import libtbx.pkg_utils

try:
    import dials.precommitbx.nagger

    dials.precommitbx.nagger.nag()
except ModuleNotFoundError:
    pass

try:
    import pkg_resources
except ModuleNotFoundError:
    pkg_resources = None


def _install_dxtbx_setup():
    """Install xia2 as a regular/editable python package"""

    # We need to find this from libtbx.env because libtbx doesn't read
    # this properly, as a module - it's just eval'd in scope
    dxtbx_root_path = libtbx.env.dist_path("dxtbx")

    # Call pip
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-build-isolation",
            "--no-deps",
            "-e",
            dxtbx_root_path,
        ],
        check=True,
    )


def _install_dxtbx_setup_readonly_fallback():
    """
    Partially install package in the libtbx build folder.

    This is a less complete installation - base python console_scripts
    entrypoints will not be installed, but the basic package metadata
    and other entrypoints should be enumerable.
    """
    dxtbx_root_path = libtbx.env.dist_path("dxtbx")
    # Find a path to "faux-install" this package to. We'll then move
    # out the .egg-link
    build_path = Path(abs(libtbx.env.build_path))
    dxtbx_build_path = build_path / "dxtbx"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--prefix",
            dxtbx_build_path,
            "--no-build-isolation",
            "--no-deps",
            "-e",
            dxtbx_root_path,
        ],
        check=True,
    )
    # Check that we can find the egg-link
    links = list(build_path.glob("**/dxtbx.egg-link"))
    if len(links) == 0:
        raise RuntimeError(
            "Could not configure dxtbx with read-only base - no .egg-link found after install"
        )
    elif len(links) > 1:
        raise RuntimeError(
            "Could not configure dxtbx with read-only base - too many .egg-link found after install"
        )

    link = links[0]
    # Move this .egg-link to the build/lib path (libtbx puts this as a root PYTHONPATH)
    dest = build_path / "lib" / link.name
    print(f"Renaming {link} → {dest}")
    link.rename(dest)


# Retain until after DIALS 3.6 is released to unregister the previous dispatcher handlers
if not pkg_resources or any(x.key == "libtbx.dxtbx" for x in pkg_resources.working_set):
    libtbx.pkg_utils.define_entry_points({})

# Try to detect case where base python environment is read-only
# e.g. on an LCLS session on a custom cctbx installation where the
# source is editable but the conda_base is read-only
#
# Pip uses this check to fallback to --user - we don't want that so check now
if os.name == "posix" and not os.access(path, os.W_OK):
    print("Python site directory not writable - falling back to tbx install")
    _install_dxtbx_setup_readonly_fallback()
else:
    # Either writable or on Windows
    _install_dxtbx_setup()
