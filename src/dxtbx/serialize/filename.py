from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import TypeVar

from dxtbx.sequence_filenames import template_string_to_glob_expr

T = TypeVar("T", str, Path)


def resolve_path(path: T, directory: str | os.PathLike | None = None) -> T:
    """Resolve a file path.

    First expand any environment and user variables. Then create the absolute
    path by applying the relative path to the provided directory, if necessary.

    Args:
      path: The path to resolve directory: The local path to resolve relative
      links

    Returns:
        The absolute path to the file to read, if accessible, otherwise
        return the original path.
    """

    assert path, "Believe this is overly defensive, check if we ever rely on"
    if not path:
        return path
    trial_path = os.path.expanduser(os.path.expandvars(str(path)))
    if directory and not os.path.isabs(trial_path):
        trial_path = os.path.join(os.fspath(directory), trial_path)
    trial_path = os.path.abspath(trial_path)
    if glob.glob(template_string_to_glob_expr(trial_path)):
        return type(path)(trial_path)
    else:
        return type(path)(path)
