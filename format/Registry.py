# A registry class to handle Format classes and provide lists of them when
# this is useful for i.e. identifying the best tool to read a given range
# of image formats.

from __future__ import absolute_import, division, print_function

import pkg_resources
import warnings


def get_format_class_for(format_class_name):
    """Return the named format class
    :param format_class_name: Name of the format class
    :return: The (uninstantiated) class object
    """
    return get_format_class_index()[format_class_name]()


def lookup(format_class):
    warnings.warn(
        "lookup() has been renamed to get_format_class_for()",
        DeprecationWarning,
        stacklevel=2,
    )
    return get_format_class_for(format_class)


def get_format_class_index():
    """Return a dictionary of all known format classes.
    :return: A dictionary containing entries
             {format class name: format class factory function}
             The factory function takes no arguments and returns an
             uninstantiated format class. This avoids importing all
             format classes.
    """
    if not hasattr(get_format_class_index, "cache"):
        setattr(
            get_format_class_index,
            "cache",
            {e.name: e.load for e in pkg_resources.iter_entry_points("dxtbx.format")},
        )
    register = get_format_class_index.cache.copy()
    return register


def get_format_class_dag():
    """Return a directed acyclical graph of the format classes.
    :return: A dictionary with entries
             {format class name: [subformat class names]}
    """
    if not hasattr(get_format_class_dag, "cache"):
        index = {name: loader() for name, loader in get_format_class_index().items()}
        dag = {}
        for name in index:
            for parent in index[name].__bases__:
                dag.setdefault(parent.__name__, []).append(name)
        for key in dag:
            dag[key].sort()
        setattr(get_format_class_dag, "cache", dag)
    dag = get_format_class_dag.cache.copy()
    return dag


_format_dag = get_format_class_dag()


def get_format_class_for_file(image_file, format_hint=None):
    """Find the best format handler in the registry for given image file
    :param image_file: A string containing the file path to an image
    :param format_hint: An optional string of a format class name that should
                        be checked first
    :return: An uninstantiated format class, if a matching one was found,
             otherwise None
    """

    # If a format hint was given then find all paths through the DAG leading
    # to this format class. Create a set containing all format class names
    # on those paths (priority_formats), and use the information whether a
    # format name is in the set to prioritise it during the format search.
    priority_formats = set()
    if format_hint:
        priority_candidates = set((format_hint,))
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
            if format_class.understand(image_file):
                return recurse(child, image_file)
        return get_format_class_for(format_name)

    # Starting at "Format" and using any potential prioritisation information
    # look for any path through the DAG of formats, stopping at the first
    # accepting leaf node
    for format in sorted(_format_dag["Format"], key=format_sort):
        format_class = get_format_class_for(format)
        if format_class.understand(image_file):
            return recurse(format, image_file)

    # There is no class accepting this file, so return None instead
    return None


class _Registry:
    # deprecated class

    @staticmethod
    def get():
        warnings.warn("Registry.get() is deprecated", DeprecationWarning, stacklevel=2)
        return tuple(
            get_format_class_for(format_class) for format_class in _format_dag["Format"]
        )

    @staticmethod
    def find(image_file, format_hint=None):
        warnings.warn("Registry.find() is deprecated", DeprecationWarning, stacklevel=2)
        return get_format_class_for_file(image_file, format_hint)


Registry = _Registry()
