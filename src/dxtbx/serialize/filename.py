import os


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
