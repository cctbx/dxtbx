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
        elif isinstance(value, list) and isinstance(value[0], list):
            for j, element in enumerate(value):
                label = "_%s/%d" % (key, j)
                if (
                    len(element) == 0
                    or isinstance(element[0], int)
                    or isinstance(element[0], float)
                ):
                    h5_ptr.create_dataset(label, data=element)
                elif isinstance(element[0], dict):
                    _element = _list_dict_to_dict_list(element)
                    _dict_to_h5(_element, h5_ptr.create_group("*%s" % label))
                elif isinstance(element[0], tuple):
                    h5_ptr.create_dataset(label, data=element)
                elif isinstance(element[0], str):
                    h5_ptr.create_dataset(label, data=numpy.array(element, dtype="S"))
                else:
                    raise TypeError(
                        "no idea what to do with [[thing]] where thing != dict, int or float: %s"
                        % str(element)
                    )

        elif isinstance(value, list) and isinstance(value[0], dict):
            _value = _list_dict_to_dict_list(value)
            _dict_to_h5(_value, h5_ptr.create_group("*%s" % key))
        elif isinstance(value, list) and isinstance(value[0], str):
            h5_ptr.create_dataset(key, data=numpy.array(value, dtype="S"))
        else:
            h5_ptr.create_dataset(key, data=value)


def dict_to_h5(data, h5_file_out, prefix):
    """Push a dictionary into a new named hdf5 file, inside a group named
    prefix, e.g. /dials/experiments."""

    with h5py.File(h5_file_out, "w") as h5_out:
        out = h5_out.create_group(prefix)
        _dict_to_h5(data, out)


def _h5_value_to_python(value):
    if hasattr(value, "shape") and len(value.shape) > 1:
        return _unpack_tuples(value)
    if type(value) == str:
        return value
    else:
        if value.dtype == numpy.dtype("int64"):
            return [int(i) for i in list(value)]
        elif value.dtype.kind == "S":
            return [s.decode() for s in value]
        else:
            return list(value)


def _unpack_tuples(value):
    if value.dtype == numpy.dtype("int64"):
        return [
            tuple([int(value[i, j]) for j in range(value.shape[1])])
            for i in range(value.shape[0])
        ]
    else:
        return [
            tuple([value[i, j] for j in range(value.shape[1])])
            for i in range(value.shape[0])
        ]


def _h5_to_dict(h5_ptr):
    result = {}
    for key in h5_ptr:
        value = h5_ptr[key]
        if isinstance(value, h5py.Group):
            if key.startswith("*"):
                result[key[1:]] = _dict_list_to_list_dict(_h5_to_dict(value))
            elif key.startswith("_"):
                tmp = []
                for v in value:
                    if isinstance(v, h5py.Group):
                        result.append(_h5_to_dict(value[v]))
                    else:
                        tmp.append(_h5_value_to_python(value[v]))
                result[key[1:]] = tmp
            else:
                result[key] = _h5_to_dict(value)
        else:
            value = value[()]
            result[key] = _h5_value_to_python(value)

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
        "list_of_lists_of_tuples": [[(j, j + 1, j + 2) for j in range(0, 12, 3)]],
    }
    dict_to_h5(things, "things.h5", "/things")
    other = h5_to_dict("things.h5", "/things")

    import json

    original = json.dumps(things, sort_keys=True)
    reloaded = json.dumps(other, sort_keys=True)

    assert original == reloaded
