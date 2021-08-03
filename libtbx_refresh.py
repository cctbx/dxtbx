import os
import random
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
    dxtbx_root_path = Path(libtbx.env.dist_path("dxtbx"))
    # Install this into a build/dxtbx subfolder
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

    # Update the libtbx environment to point to the source location
    libtbx.env.pythonpath.insert(0, str(dxtbx_root_path / "src"))
    libtbx.env.pickle()


def _test_writable_dir(path: Path) -> bool:
    """Test a path is writable. Based on pip's _test_writable_dir_win."""
    # os.access doesn't work on windows
    # os.access won't always work with network filesystems
    # pip doesn't use tempfile on windows because https://bugs.python.org/issue22107
    basename = "test_site_packages_writable_dxtbx_"
    alphabet = "abcdefghijklmnopqrstuvwxyz0123456789"
    for _ in range(10):
        name = basename + "".join(random.choice(alphabet) for _ in range(6))
        file = path / name
        try:
            fd = os.open(file, os.O_RDWR | os.O_CREAT | os.O_EXCL)
        except FileExistsError:
            pass
        except PermissionError:
            return False
        else:
            os.close(fd)
            os.unlink(file)
            return True


site_packages = Path(site.getsitepackages()[0])

# Detect case where base python environment is read-only
# e.g. on an LCLS session on a custom cctbx installation where the
# source is editable but the conda_base is read-only
#
# We need to check before trying to install as pip does os.access-based
# checks then installs with --user if it fails. We don't want that.
if _test_writable_dir(site_packages):
    _install_dxtbx_setup()
else:
    print("Python site directory not writable - falling back to tbx install")
    _install_dxtbx_setup_readonly_fallback()

# Retain until after DIALS 3.6 is released to unregister the previous dispatcher handlers
if not pkg_resources or any(x.key == "libtbx.dxtbx" for x in pkg_resources.working_set):
    libtbx.pkg_utils.define_entry_points({})
