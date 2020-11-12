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

plugin_path = libtbx.env.under_base(os.path.join("lib", "plugins"))
new_plugin_path = libtbx.env.under_base(os.path.join("lib", "hdf5", "plugin"))
if os.path.exists(plugin_path) and not os.path.exists(new_plugin_path):
    # Set up the plugin path for HDF5 to pick up compression plugins from legacy location
    os.environ["HDF5_PLUGIN_PATH"] = (
        plugin_path + os.pathsep + os.getenv("HDF5_PLUGIN_PATH", "")
    )
    warnings.warn(
        "You are using an outdated version of the hdf5-external-filter-plugins package.",
        DeprecationWarning,
    )

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
