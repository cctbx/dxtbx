from __future__ import absolute_import, division, print_function

import os
import warnings


def load_path(path, directory=None):
    """[DEPRECATED: Use resolve_path] Load a filename from a JSON file.

    First expand any environment and user variables. Then create the absolute path
    from the current directory (i.e. the directory in which the JSON file is
    located.

    Params:
      path The path to the file.

    """
    warnings.warn(
        "Use dxtbx.filename.resolve_path instead of load_path",
        DeprecationWarning,
        stacklevel=2,
    )
    return resolve_path(path, directory)


def resolve_path(path, directory=None):
    """Resolve a file path.

    First expand any environment and user variables. Then create the absolute
    path by applying the relative path to the provided directory, if necessary.

    Args:
      path (str): The path to resolve
      directory (Optional[str]): The local path to resolve relative links

    Returns:
        str: The absolute path to the file to read
    """
    if not path:
        return ""
    path = os.path.expanduser(os.path.expandvars(path))
    if directory and not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)
