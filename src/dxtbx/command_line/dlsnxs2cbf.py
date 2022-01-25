"""
Convert a NXmx-format NeXus file to a set of CBF-format image files.

Note that this tool does not produce full imgCIF-format files, only
Dectris-style mini-CBF files consisting of a plain text simplified
header and the binary compressed image data.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import dxtbx.util.dlsnxs2cbf

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nexus_file", metavar="nexus-file", help="Input NeXus file.", type=Path
)
parser.add_argument(
    "-o",
    "--output-directory",
    help="Directory in which to store the CBF image files.  Defaults to the current "
    "working directory.",
    type=Path,
    default=Path(),
)
parser.add_argument(
    "-p",
    "--parents",
    help="Create all parents of the output directory if they do not already exist.  "
    "By default, an exception will be raised if part of the output path is missing.",
    action="store_true",
)
parser.add_argument(
    "-t",
    "--name-template",
    help="Template for the CBF file name stem.  The output filenames will be "
    "constructed by appending an image number and '.cbf' file extension to this "
    "template.  For example, the template 'image_' will result in files named like "
    "'image_001.cbf'.  To name the files by image number alone, e.g. '0001.cbf', "
    "pass an empty template string with '-t \"\"'.  The default template is the stem "
    "of the input NeXus file name, with an added trailing underscore.",
)
parser.add_argument(
    "-n",
    "--number-of-digits",
    help="The number of digits to use for the image number label in the CBF file "
    "names.  If this exceeds the minimal necessary number of digits to accommodate "
    "the number of the last CBF image, the image number label will be padded with "
    "zeros.  If an insufficient number is specified, the minimal necessary number of "
    "digits will be used instead, which is the default behaviour.",
    type=int,
    default=0,
)


def run(args: list[str] | None = None):
    dxtbx.util.encode_output_as_utf8()

    args = parser.parse_args(args)

    try:
        args.output_directory.mkdir(parents=args.parents, exist_ok=True)
    except FileNotFoundError as e:
        sys.exit(e)

    dxtbx.util.dlsnxs2cbf.make_cbf(
        args.nexus_file,
        args.output_directory,
        args.name_template,
        args.number_of_digits,
    )


if __name__ == "__main__":
    run()
