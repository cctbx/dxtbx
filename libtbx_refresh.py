import subprocess
import sys

import libtbx.load_env

try:
    import dials.precommitbx.nagger

    dials.precommitbx.nagger.nag()
except ModuleNotFoundError:
    pass


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


_install_dxtbx_setup()
