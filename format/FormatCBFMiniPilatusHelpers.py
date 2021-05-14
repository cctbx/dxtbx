"""Helpers for FormatCBFMiniPilatus..."""


import calendar
import time


def get_pilatus_timestamp(timestamp_string):
    if "." in timestamp_string:
        timestamp, milliseconds = timestamp_string.split(".")
    else:
        timestamp = timestamp_string
        milliseconds = "000"

    for format in ["%Y-%b-%dT%H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y/%b/%d %H:%M:%S"]:

        try:
            struct_time = time.strptime(timestamp, format)
            return calendar.timegm(struct_time) + float("0." + milliseconds)

        except Exception:
            pass

    raise RuntimeError("timestamp %s not recognised" % timestamp)
