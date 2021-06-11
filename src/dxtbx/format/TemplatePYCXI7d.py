import sys

from dxtbx.format.FormatPYunspecified import FormatPYunspecified


class FormatPYCXI7d(FormatPYunspecified):

    """PREFERENCE FILE.
    Treats any Pickle-format file lacking a DETECTOR_FORMAT_VERSION key
    automatically as format CXI 7.d.
    Rename this file to FormatPYCXI7d.py and put in dxtbx/format path.
    """

    def _start(self):

        FormatPYunspecified.start_helper(
            self, version_token="distl.detector_format_version=CXI 7.d"
        )


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(FormatPYCXI7d.understand(arg))
