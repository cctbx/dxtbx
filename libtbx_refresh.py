from __future__ import annotations

import contextlib
import importlib
import inspect
import io
import os
import subprocess
import sys
from pathlib import Path

import libtbx
import libtbx.pkg_utils

try:
    import dials.precommitbx.nagger

    dials.precommitbx.nagger.nag()
except ModuleNotFoundError:
    pass


def _find_site_packages_with_metadata(package_name: str, build_path: Path):
    """
    Find the site-packages directory containing the package metadata.
    Returns the site-packages directory if metadata is found, None otherwise.
    """
    # Look for Python site-packages directories in the build path
    for python_dir in build_path.glob("lib/python*"):
        site_packages = python_dir / "site-packages"
        if site_packages.exists():
            # Look for both .dist-info and .egg-info directories
            for pattern in [f"{package_name}*.dist-info", f"{package_name}*.egg-info"]:
                for metadata_dir in site_packages.glob(pattern):
                    if metadata_dir.exists():
                        # Return site-packages only if we actually found metadata
                        return site_packages
            # If no metadata found in this site-packages, continue searching
    return None


def _install_setup_readonly_fallback(package_name: str):
    """
    Partially install package in the libtbx build folder.

    This is a less complete installation - base python console_scripts
    entrypoints will not be installed, but the basic package metadata
    and other entrypoints will be enumerable through dispatcher black magic
    """
    root_path = libtbx.env.dist_path(package_name)
    import_path = os.path.join(root_path, "src")

    # Install this into a build/dxtbx subfolder
    build_path = abs(libtbx.env.build_path / package_name)
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--prefix",
            build_path,
            "--no-build-isolation",
            "--no-deps",
            "-e",
            root_path,
        ],
        check=True,
    )

    # Get the actual environment being configured (NOT libtbx.env)
    env = _get_real_env_hack_hack_hack()

    # As of PEP 660, the package metadata (dist-info) goes in the install dir,
    # not the source dir. Add this location to the python path.
    metadata_dir = _find_site_packages_with_metadata(package_name, Path(build_path))
    if metadata_dir and metadata_dir not in sys.path:
        metadata_rel = libtbx.env.as_relocatable_path(str(metadata_dir))
        if metadata_rel not in env.pythonpath:
            env.pythonpath.insert(0, metadata_rel)

    # Update the sys.path so we can find the package in this process
    # if we do a full reconstruction of the working set
    if import_path not in sys.path:
        sys.path.insert(0, import_path)

    # ...and make sure it is picked up by the import system
    importlib.invalidate_caches()

    # This is already generated by this point, but will get picked up
    # on the second libtbx.refresh.
    module = env.module_dict[package_name]
    if f"src/{package_name}" not in module.extra_command_line_locations:
        module.extra_command_line_locations.append(f"src/{package_name}")

    # Because dispatchers for all modules are generated _before_ any of
    # libtbx_refresh are run, then we need to regenerate all of the
    # dispatchers now we've added the extra PYTHONPATH
    with contextlib.redirect_stdout(io.StringIO()):
        for module in env.module_list:
            module.process_command_line_directories()


def _get_real_env_hack_hack_hack():
    """
    Get the real, currently-being-configured libtbx.env environment.

    This is not libtbx.env, because although libtbx.env_config.environment.cold_start
    does:
        self.pickle()
        libtbx.env = self
    the first time there is an "import libtbx.load_env" this environment
    gets replaced by unpickling the freshly-written libtbx_env file onto
    libtbx.env, thereby making the environment accessed via libtbx.env
    *not* the actual one that is currently being constructed.

    So, the only way to get this environment being configured in order
    to - like - configure it, is to walk the stack trace and extract the
    self object from environment.refresh directly.
    """
    for frame in inspect.stack():
        if (
            frame.filename.endswith("env_config.py")
            and frame.function == "refresh"
            and "self" in frame.frame.f_locals
        ):
            return frame.frame.f_locals["self"]

    raise RuntimeError("Could not determine real libtbx.env_config.environment object")


# When building in libtbx, always assume it's unsafe to write to base/
_install_setup_readonly_fallback("dxtbx")
