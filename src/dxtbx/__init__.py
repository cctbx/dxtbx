from __future__ import annotations

import logging
import os
import sys
import typing

import libtbx.load_env

import dxtbx.format.Registry

if sys.version_info.major == 2:
    sys.exit("Python 2 is no longer supported")

from .version import version as __version__  # noqa: F401

if typing.TYPE_CHECKING:
    from dxtbx.format.Format import Format

# Ensures that HDF5 has the conda_base plugin path configured.
#
# Ideally this will be properly configured by the conda environment.
# However, currently the dials-installer will not install a path-correct
# conda_base folder, so it needs to be updated manually.

_hdf5_plugin_path = libtbx.env.under_base(os.path.join("lib", "hdf5", "plugin"))

# Inject via the environment if h5py not used yet, or else use h5py
if "h5py" not in sys.modules:
    os.environ["HDF5_PLUGIN_PATH"] = (
        _hdf5_plugin_path + os.pathsep + os.getenv("HDF5_PLUGIN_PATH", "")
    )
else:
    # We've already loaded h5py, so setting the environment variable won't work
    # We need to use the h5py API to add a plugin path
    import h5py

    h5_plugin_paths = [h5py.h5pl.get(i).decode() for i in range(h5py.h5pl.size())]
    if _hdf5_plugin_path not in h5_plugin_paths:
        h5py.h5pl.prepend(_hdf5_plugin_path.encode())

logging.getLogger("dxtbx").addHandler(logging.NullHandler())


class IncorrectFormatError(RuntimeError):
    """
    An exception class for an incorrect format error
    """

    def __init__(self, format_instance, filename):
        super().__init__(
            "Could not open {} as {}".format(filename, str(format_instance))
        )
        self.args = (format_instance, filename)


def load(filename: str | bytes | os.PathLike) -> Format:
    """Use DXTBX to load the input filename.

    Args:
        filename:  The input filename

    Returns
        A dxtbx Format-subclass instance for the file type

    Raises:
        IOError: If the file format could not be determined
    """
    format_instance = dxtbx.format.Registry.get_format_class_for_file(filename)
    return format_instance(os.fspath(filename))
