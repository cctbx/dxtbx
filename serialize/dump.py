from __future__ import absolute_import, division, print_function

import json
import re
import textwrap
import warnings


def compact_simple_list(match):
    """Callback function. Given a simple list match, compact it and ensure
    that it wraps around by 80 characters.

    Params:
        match The regular expression match

    Returns:
        The string to replace the expression with

    """
    # Calculate the initial indent as the length of the first match group
    initial_indent = len(match.group(1))

    # Get the lines in the match
    lines = match.group(2).splitlines()

    # Set the indent by finding the indent of the first lines
    if len(lines) > 1:
        subsequent_indent = len(lines[1]) - len(lines[1].lstrip())
    else:
        subsequent_indent = 0

    # Strip whitespace from the lines
    lines = [l.strip() for l in lines]

    # Create and return the string wrapped about 80 chars
    list_string = "\n".join(
        textwrap.wrap(
            " ".join(lines),
            80,
            initial_indent=" " * initial_indent,
            subsequent_indent=" " * subsequent_indent,
        )
    ).lstrip()

    # Return the string
    return match.group(1) + list_string


def compact_simple_lists(string):
    """Find simple lists in the string and compact.

    Params:
        string The input JSON string

    Returns:
        The output JSON string

    """
    return re.sub(r'(.*"\w+".*:.*)(\[[^\{\}\[\]]*\])', compact_simple_list, string)


def imageset_to_string(obj, compact=False):
    """Dump the given object to string.

    Params:
        obj The imageset
        compact Write in compact representation

    Returns:
        The JSON string

    """
    from dxtbx.serialize.imageset import imageset_to_dict

    if compact:
        return json.dumps(
            imageset_to_dict(obj), separators=(",", ":"), ensure_ascii=True
        )
    else:
        return json.dumps(imageset_to_dict(obj), indent=2, ensure_ascii=True)


def imageset(obj, outfile, compact=False):
    """Dump the given object to file.

    Params:
        obj The imageset to dump
        outfile The output file name or file object
        compact Write in compact representation

    """
    # If the input is a string then open and write to that file
    if isinstance(outfile, str):
        with open(outfile, "wb") as outfile:
            outfile.write(imageset_to_string(obj, compact).encode())

    # Otherwise assume the input is a file and write to it
    else:
        outfile.write(imageset_to_string(obj, compact).encode())


def datablock(obj, outfile, **kwargs):
    """ Dump the given object to file. """
    from dxtbx.datablock import DataBlockDumper

    dump = DataBlockDumper(obj)
    dump.as_file(outfile, **kwargs)


def experiment_list(obj, outfile):
    """ Dump an experiment list. """
    warnings.warn(
        "use .as_file() on the experimentlist directly",
        DeprecationWarning,
        stacklevel=2,
    )

    obj.as_file(outfile)
