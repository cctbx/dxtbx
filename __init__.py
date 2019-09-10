from __future__ import absolute_import, division, print_function

import logging
import os
import sys
import warnings

import dxtbx.format.Registry

# Invert FPE trap defaults, https://github.com/cctbx/cctbx_project/pull/324
if "boost.python" in sys.modules:
    import boost.python

    boost.python.ext.trap_exceptions(
        bool(os.getenv("BOOST_ADAPTBX_TRAP_FPE")),
        bool(os.getenv("BOOST_ADAPTBX_TRAP_INVALID")),
        bool(os.getenv("BOOST_ADAPTBX_TRAP_OVERFLOW")),
    )
elif not os.getenv("BOOST_ADAPTBX_TRAP_FPE") and not os.getenv(
    "BOOST_ADAPTBX_TRAP_OVERFLOW"
):
    os.environ["BOOST_ADAPTBX_FPE_DEFAULT"] = "1"


def _deprecate_function(name, func):
    def _wrapper(*args, **kwargs):
        warnings.warn(
            "Addressing dxtbx extensions as dxtbx.{0} is deprecated. Instead use dxtbx.ext.{0}".format(
                name
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        return func(*args, **kwargs)

    return _wrapper


try:
    import dxtbx.ext as _ext

    for funcname in dir(_ext):
        if funcname != "ext" and not funcname.startswith("_"):
            globals()[funcname] = _deprecate_function(funcname, getattr(_ext, funcname))
except ImportError:
    pass

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
