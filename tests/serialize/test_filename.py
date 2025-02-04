from __future__ import annotations

import os
import uuid

from dxtbx.serialize.filename import resolve_path


def test_resolve_path(monkeypatch):
    # Resolve path
    def alwaystrue(path):
        return True

    # Set an environment variable
    monkeypatch.setenv("HELLO_WORLD", "EXPANDED")

    # First try a path that does not exist
    filename = str(uuid.uuid4())
    new_path = os.path.join("~", "$HELLO_WORLD", filename)
    path = resolve_path(new_path)
    assert path == new_path

    # Now pretend the path exists
    monkeypatch.setattr(os.path, "exists", alwaystrue)
    path = resolve_path(new_path)
    assert path == os.path.join(os.path.expanduser("~"), "EXPANDED", filename)

    new_path = os.path.join("$HELLO_WORLD", filename)
    path = resolve_path(new_path)
    assert path == os.path.abspath(os.path.join("EXPANDED", filename))
