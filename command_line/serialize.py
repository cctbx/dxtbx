from __future__ import absolute_import, division, print_function

from optparse import OptionParser

from dxtbx.imageset import ImageSetFactory
from dxtbx.serialize import dump

if __name__ == "__main__":
    # Specify the command line options
    usage = "usage: %prog [options] /path/to/image/files.ext"
    parser = OptionParser(usage)

    # Add a verbose option (False by default)
    parser.add_option(
        "-o",
        "--output-file",
        dest="output_file",
        type="string",
        default="imageset.json",
        help="Enter a destination filename for serialization",
    )

    # Parse the arguments
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) == 0:
        print(parser.print_help())

    else:
        imagesets = ImageSetFactory.new(args)
        if len(imagesets) == 0:
            print("Error: no imagesets to serialize.")
        elif len(imagesets) > 1:
            print("Error: more than 1 imageset has been specified")
        else:
            dump.imageset(imagesets[0], options.output_file)
            print("Serialized imageset to {}".format(options.output_file))
