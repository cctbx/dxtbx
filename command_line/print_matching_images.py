import argparse

from dxtbx.sequence_filenames import find_matching_images


def run(args=None):
    parser = argparse.ArgumentParser(
        description="Find images that match a template specification"
    )
    parser.add_argument(
        "template", help="The template specification", metavar="TEMPLATE"
    )
    options = parser.parse_args(args)

    for mi in find_matching_images(options.template):
        print(mi)


if __name__ == "__main__":
    run()
