from __future__ import division, print_function

import h5py
import numpy

__all__ = ["dict_to_h5", "h5_to_dict"]

# cannot store compound data types e.g. list of dicts so convert this to a
# dict of lists which we can store. In application the need for this inversion
# will be tagged by prefixing the dataset name with * which seems a Pythonic
# convention


def _list_dict_to_dict_list(list_of_dict):
    return {k: [d[k] for d in list_of_dict] for k in list_of_dict[0]}


def _dict_list_to_list_dict(dict_of_list):
    return [dict(zip(dict_of_list, t)) for t in zip(*dict_of_list.values())]


def _dict_to_h5(data, h5_ptr):
    for key in sorted(data):
        value = data[key]
        if isinstance(value, dict):
            _dict_to_h5(value, h5_ptr.create_group(key))
        elif isinstance(value, list) and len(value) == 0:
            h5_ptr.create_dataset(key, data=value)
        elif isinstance(value, list) and isinstance(value[0], dict):
            _value = _list_dict_to_dict_list(value)
            _dict_to_h5(_value, h5_ptr.create_group("*%s" % key))
        elif isinstance(value, list) and isinstance(value[0], str):
            h5_ptr.create_dataset(key, data=numpy.array(value, dtype="S"))
        else:
            print(key, value)
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
            if key.startswith("*"):
                result[key[1:]] = _dict_list_to_list_dict(_h5_to_dict(value))
            else:
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
        "list_of_dicts": [{"name": "foo", "value": j} for j in range(10)],
        "list_of_lists": [list(range(j)) for j in range(5, 10)],
    }
    dict_to_h5(things, "things.h5", "/things")
    other = h5_to_dict("things.h5", "/things")
