from __future__ import absolute_import, division, print_function

import math

from future.builtins import round


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
        fmt_str = "%%.%if(%%i)" % precision
        return fmt_str % (value, round(su))
    else:
        precision += 1
        su = int(round(standard_uncertainty, precision))
        fmt_str = "%.0f(%i)"
        return fmt_str % (round(value, precision), su)
