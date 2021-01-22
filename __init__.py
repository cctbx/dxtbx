from __future__ import absolute_import, division, print_function

import logging
import os
import sys
import warnings

import libtbx.load_env

import dxtbx.format.Registry

if sys.version_info.major == 2:
    warnings.warn(
        "Python 2 is no longer supported. "
        "If you need Python 2.7 support please use the DIALS 2.2 release branch.",
        UserWarning,
    )


# Ensures that HDF5 has the conda_base plugin path configured.
#
# Ideally this will be properly configured by the conda environment.
# However, currently the dials-installer will not install a path-correct
# conda_base folder, so it needs to be updated manually.

# Warn if user is using an old hdf5-external-filter-plugins
# Remove after DIALS 3.3 release branch is made
_legacy_plugin_path = libtbx.env.under_base(os.path.join("lib", "plugins"))
_new_plugin_path = libtbx.env.under_base(os.path.join("lib", "hdf5", "plugin"))

# Find out which plugin path we need to inject
if os.path.exists(_legacy_plugin_path) and not os.path.exists(_new_plugin_path):
    # Set up the plugin path for HDF5 to pick up compression plugins from legacy location
    _hdf5_plugin_path = _legacy_plugin_path
    warnings.warn(
        "You are using an outdated version of the hdf5-external-filter-plugins package.\n"
        "Please update your environment using\n"
        "'conda install \"hdf5-external-filter-plugins==0.1.0[build_number='>=5']\"'",
        UserWarning,
    )
else:
    _hdf5_plugin_path = _new_plugin_path

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
        super(IncorrectFormatError, self).__init__(
            "Could not open %s as %s" % (filename, str(format_instance))
        )
        self.args = (format_instance, filename)


def load(filename):
    """Use DXTBX to load the input filename.

    :param filename:  The input filename
    :type  filename:  os.PathLike or str or bytes
    :returns:         A dxtbx Format-subclass instance for the file type
    :raises IOError:  if the file format could not be determined
    """
    # Unwrap PEP-519-style objects. This covers py.path, pathlib, ...
    if hasattr(filename, "__fspath__"):
        filename = filename.__fspath__()
    format_instance = dxtbx.format.Registry.get_format_class_for_file(filename)
    return format_instance(filename)
