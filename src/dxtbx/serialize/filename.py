from __future__ import annotations

import glob
import os
from typing import AnyStr

from dxtbx.sequence_filenames import template_string_to_glob_expr


def resolve_path(
    path: AnyStr | os.PathLike, directory: AnyStr | os.PathLike | None = None
) -> str:
    """Resolve a file path.

    First expand any environment and user variables. Then create the absolute
    path by applying the relative path to the provided directory, if necessary.

    Args:
      path (str): The path to resolve
      directory (Optional[str]): The local path to resolve relative links

    Returns:
        str: The absolute path to the file to read if accessible, otherwise
        return the original path as provided
    """

    path = str(path)
    if not path:
        return ""
    trial_path = os.path.expanduser(os.path.expandvars(path))
    if directory and not os.path.isabs(trial_path):
        trial_path = os.path.join(str(directory), trial_path)
    trial_path = os.path.abspath(trial_path)
    if glob.glob(template_string_to_glob_expr(trial_path)):
        return trial_path
    else:
        return path
