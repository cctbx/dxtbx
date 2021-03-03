import argparse

import dxtbx.util.dlsnxs2cbf


def run(args=None):
    dxtbx.util.encode_output_as_utf8()
    parser = argparse.ArgumentParser(description="Convert an nxs file to multiple cbfs")
    parser.add_argument("nexus_file", help="Input nexus file")
    parser.add_argument(
        "template", help="Template cbf output name e.g. 'image_%04d.cbf'"
    )
    options = parser.parse_args(args)
    dxtbx.util.dlsnxs2cbf.make_cbf(options.nexus_file, options.template)


if __name__ == "__main__":
    run()
