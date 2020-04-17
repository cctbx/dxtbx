from __future__ import division, print_function

import h5py
import numpy

__all__ = ["dict_to_h5", "h5_to_dict"]


def _dict_to_h5(data, h5_ptr):
    for key in sorted(data):
        value = data[key]
        if isinstance(value, dict):
            _dict_to_h5(value, h5_ptr.create_group(key))
        else:
            h5_ptr.create_dataset(key, data=value)


def dict_to_h5(data, h5_file_out, prefix):
    """Push a dictionary into a new named hdf5 file, inside a group named
    prefix, e.g. /dials/experiments."""

    with h5py.File(h5_file_out, "w") as h5_out:
        out = h5_out.create_group(prefix)
        _dict_to_h5(data, out)


def _h5_to_dict(h5_ptr):
    result = {}
    for key in h5_ptr:
        value = h5_ptr[key]
        if isinstance(value, h5py.Group):
            result[key] = _h5_to_dict(value)
        else:
            value = value[()]
            if type(value) == str:
                result[key] = value
            else:
                if value.dtype == numpy.dtype("int64"):
                    result[key] = [int(i) for i in list(value)]
                else:
                    result[key] = list(value)
    return result


def h5_to_dict(h5_file_in, prefix):
    """Read the dataset named by prefix and extract the contents as a dictionary
    - do not do this with e.g. an Eiger dataset unless you have a lot of
    time and memory."""

    with h5py.File(h5_file_in, "r") as h5_in:
        return _h5_to_dict(h5_in[prefix])


if __name__ == "__main__":
    import string

    things = {
        "ints": list(range(100)),
        "floats": list(0.01 * j for j in range(100)),
        "letters": {"upper": string.ascii_uppercase, "lower": string.ascii_lowercase},
        "keys": dict(zip(string.ascii_uppercase, string.ascii_lowercase)),
    }
    dict_to_h5(things, "things.h5", "/things")
    other = h5_to_dict("things.h5", "/things")

    # compare a dumb way

    import json

    original = json.dumps(things, sort_keys=True)
    reloaded = json.dumps(other, sort_keys=True)

    assert original == reloaded
