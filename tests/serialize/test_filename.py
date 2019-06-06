from __future__ import absolute_import, division, print_function

import pytest
import os

from dxtbx.serialize.filename import load_path, resolve_path


def test_resolve_path():
    os.environ["HELLO_WORLD"] = "EXPANDED"
    new_path = os.path.join("~", "$HELLO_WORLD", "path")
    path = resolve_path(new_path)
    assert path == os.path.join(os.path.expanduser("~"), "EXPANDED", "path")
    new_path = os.path.join("$HELLO_WORLD", "path")
    path = resolve_path(new_path)
    assert path == os.path.abspath(os.path.join("EXPANDED", "path"))


def test_load_path_deprecated():
    os.environ["HELLO_WORLD"] = "EXPANDED"
    new_path = os.path.join("~", "$HELLO_WORLD", "path")
    with pytest.deprecated_call():
        path = load_path(new_path)
    assert path == resolve_path(new_path)
