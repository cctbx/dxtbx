from __future__ import annotations

import json
import math
import os
import sys
from urllib.parse import urlparse

import libtbx


class AutoEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, libtbx.AutoType):
            return "Auto"
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


def encode_output_as_utf8() -> None:
    """
    Ensures stdout and stderr are UTF-8 encoded.
    This avoids crashes on Windows when printing unicode where the output is
    redirected away from the console into a pipe or file
    """
    if hasattr(encode_output_as_utf8, "done"):
        return

    if sys.stdout.encoding.lower() != "utf-8":
        sys.stdout.reconfigure(encoding="utf-8")  # type: ignore
        sys.stderr.reconfigure(encoding="utf-8")  # type: ignore
    encode_output_as_utf8.done = True  # type: ignore


def format_float_with_standard_uncertainty(value, standard_uncertainty, minimum=1e-12):
    """
    Formats a float, including the uncertainty in its value.

    Parameters
    ----------
    value : float
    standard_uncertainty : float
    minimum : float

    Returns
    -------
    str

    Examples
    --------
    >>> format_float_with_standard_uncertainty(5e-3, 1e-3)
    '0.0050(10)'
    >>> format_float_with_standard_uncertainty(5e-3, 1e-6)
    '0.0050000(10)'
    """
    if standard_uncertainty <= minimum:
        dp = -int(math.log10(minimum))
        return str(round(value, dp))
    precision = -int(round(math.log10(standard_uncertainty)))
    if precision > -1:
        su = standard_uncertainty * math.pow(10, precision)
        if round(su, 1) < 2:
            su *= 10
            precision += 1
        return "{value:.{precision}f}({irsu})".format(
            value=value, precision=precision, irsu=int(round(su))
        )
    else:
        precision += 1
        su = int(round(standard_uncertainty, precision))
        fmt_str = "%.0f(%i)"
        return fmt_str % (round(value, precision), su)


def show_mask_info(expt_list):
    """Print the number of masked pixels for each module in the detector,
    for each experiment in the input list."""
    for i in expt_list.imagesets():
        d = i.get_detector()
        m = i.get_mask(0)
        print("---- ----")
        print(d)
        for j, _m in enumerate(m):
            print(
                "Module {} has {} masked pixels of {}".format(
                    j, _m.count(False), _m.size()
                )
            )


def get_url_scheme(url) -> str:
    """Extract the URL scheme from the string url, respecting Windows file paths"""

    return urlparse(str(url)).scheme if "://" in str(url) else ""


def ersatz_uuid4() -> str:
    """Generate an ersatz UUID4 i.e. random 128 bit UUID."""

    h = "%032x" % int.from_bytes(os.urandom(16), byteorder="little")
    return "%s-%s-%s-%s-%s" % (h[:8], h[8:12], h[12:16], h[16:20], h[20:])
