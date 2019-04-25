from __future__ import absolute_import, division, print_function

import os

from dxtbx.serialize.filename import resolve_path


def test_resolve_path():
    os.environ["HELLO_WORLD"] = "EXPANDED"
    new_path = os.path.join("~", "$HELLO_WORLD", "path")
    path = resolve_path(new_path)
    assert path == os.path.join(os.path.expanduser("~"), "EXPANDED", "path")
    new_path = os.path.join("$HELLO_WORLD", "path")
    path = resolve_path(new_path)
    assert path == os.path.abspath(os.path.join("EXPANDED", "path"))
