import sys

from dxtbx.format.FormatPYunspecified import FormatPYunspecified


class FormatPYCXI71(FormatPYunspecified):

    """PREFERENCE FILE.
    Treats any Pickle-format file lacking a DETECTOR_FORMAT_VERSION key
    automatically as format CXI 7.1.
    Rename this file to FormatPYCXI71.py and put in dxtbx/format path.
    """

    def _start(self):

        FormatPYunspecified.start_helper(
            self, version_token="distl.detector_format_version=CXI 7.1"
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPYCXI71.understand(arg))
