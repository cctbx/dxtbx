import subprocess
import sys

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
        [sys.executable, "-m", "pip", "install", "--no-deps", "-e", dxtbx_root_path],
        check=True,
    )


# Retain until after DIALS 3.6 is release to unregister the previous dispatcher handlers
if not pkg_resources or any(x.key == "libtbx.dxtbx" for x in pkg_resources.working_set):
    libtbx.pkg_utils.define_entry_points({})

_install_dxtbx_setup()
