"""
A registry class to handle Format classes and provide lists of them when
this is useful for i.e. identifying the best tool to read a given range
of image formats.
"""
from __future__ import annotations

import os
import typing
from typing import Callable

import pkg_resources

from dxtbx.util import get_url_scheme

if typing.TYPE_CHECKING:
    from dxtbx.format.Format import Format


def get_format_class_for(format_class_name: str) -> type[Format]:
    """Return the named format class
    :param format_class_name: Name of the format class
    :return: The (uninstantiated) class object
    """
    return get_format_class_index()[format_class_name][0]()


def get_format_class_index() -> dict[str, tuple[Callable[[], type[Format]], list[str]]]:
    """Return a dictionary of all known format classes.
    :return: A dictionary containing entries
             {format_class_name: (format_class_factory_function, [base_class_names])}
             The factory function takes no arguments and returns an
             uninstantiated format class. This avoids importing all
             format classes.
    """
    if not hasattr(get_format_class_index, "cache"):
        class_index = {}
        for e in pkg_resources.iter_entry_points("dxtbx.format"):
            if ":" in e.name:
                format_name, base_classes_str = e.name.split(":", 1)
                base_classes = tuple(base_classes_str.split(","))
            else:
                format_name, base_classes = e.name, ()
            class_index[format_name] = (e.load, base_classes)
        setattr(get_format_class_index, "cache", class_index)
    register = get_format_class_index.cache.copy()  # type: ignore
    return register


def get_format_class_dag() -> dict[str, list[str]]:
    """Return a directed acyclical graph of the format classes.
    :return: A dictionary with entries
             {format class name: [subformat class names]}
    """
    if not hasattr(get_format_class_dag, "cache"):
        index = {
            name: class_info[1] for name, class_info in get_format_class_index().items()
        }
        dag: dict[str, list[str]] = {}
        for name in index:
            for parent in index[name]:
                dag.setdefault(parent, []).append(name)
        for key in dag:
            dag[key].sort()
        setattr(get_format_class_dag, "cache", dag)
    dag = get_format_class_dag.cache.copy()  # type: ignore
    return dag


_format_dag: dict[str, list[str]] = get_format_class_dag()


def get_format_class_for_file(
    image_file: str | bytes | os.PathLike, format_hint: str | None = None
) -> type[Format] | None:
    """Find the best format handler in the registry for given image file
    :param image_file: A string containing the file path to an image
    :param format_hint: An optional string of a format class name that should
                        be checked first
    :return: An uninstantiated format class, if a matching one was found,
             otherwise None
    """

    # Grab the scheme from this URI
    scheme = get_url_scheme(image_file)

    image_file_str = os.fspath(image_file)

    # If a format hint was given then find all paths through the DAG leading
    # to this format class. Create a set containing all format class names
    # on those paths (priority_formats), and use the information whether a
    # format name is in the set to prioritise it during the format search.
    priority_formats = set()
    if format_hint:
        priority_candidates = {format_hint}
        dagset = {f: set(subf) for f, subf in _format_dag.items()}
        while priority_candidates:
            priority_formats.update(priority_candidates)
            priority_candidates = {
                f
                for f, subf in dagset.items()
                if priority_candidates.intersection(subf)
            }

    def format_sort(format_name):
        return (0 if format_name in priority_formats else 1, format_name)

    def recurse(format_name, image_file):
        # Recursively check whether any of the children understand image_file,
        # in which case they are preferred over the parent format.
        for child in sorted(_format_dag.get(format_name, []), key=format_sort):
            format_class = get_format_class_for(child)
            if scheme in format_class.schemes and format_class.understand(image_file):
                return recurse(child, image_file)
        return get_format_class_for(format_name)

    # Starting at "Format" and using any potential prioritisation information
    # look for any path through the DAG of formats, stopping at the first
    # accepting leaf node
    for format in sorted(_format_dag["Format"], key=format_sort):
        format_class = get_format_class_for(format)
        if scheme in format_class.schemes and format_class.understand(image_file_str):
            return recurse(format, image_file_str)

    return None
