"""
Testing URI-passing in the Format class hierarchy
"""

from typing import Type

import pytest

import dxtbx
import dxtbx.format.Registry as Registry
from dxtbx.format.Format import Format


@pytest.fixture
def registry(monkeypatch):
    """Temporarily register a format class"""
    _temporary_formats = {}
    _original_get_format_class_for = Registry.get_format_class_for

    def _get_format_class_for(format_class_name: str) -> Type[Format]:
        """Intercept fetching format classes"""
        if format_class_name in _temporary_formats:
            return _temporary_formats[format_class_name]
        else:
            return _original_get_format_class_for(format_class_name)

    monkeypatch.setattr(Registry, "get_format_class_for", _get_format_class_for)

    class _registry:
        @staticmethod
        def register(format_class: Format) -> None:
            name = format_class.__name__
            _temporary_formats[name] = format_class
            monkeypatch.setitem(Registry._format_dag, name, [])
            _format_entry = list(Registry._format_dag["Format"]) + [name]
            monkeypatch.setitem(Registry._format_dag, "Format", _format_entry)
            # Registry._format_dag["Format"].append(name)

    return _registry


def test_registry_inject(registry, tmp_path):
    """meta-test of registry injection"""
    _hit_understand = []

    class TestClass(Format):
        @classmethod
        def understand(cls, filename):
            _hit_understand.append(True)
            return "TEST_FILE" in str(filename)

    filename_notest = tmp_path / "A_FILE"
    filename_notest.touch()
    filename = tmp_path / "TEST_FILE"
    filename.touch()
    assert not Registry.get_format_class_for_file(str(filename_notest))
    registry.register(TestClass)
    assert Registry.get_format_class_for("TestClass") is TestClass
    assert Registry.get_format_class_for_file(str(filename))
    assert _hit_understand == [True]


def test_custom_scheme_handler(registry):
    class SchemeHandler(Format):
        schemes = ["scheme"]

        @classmethod
        def understand(cls, endpoint):
            assert endpoint.startswith("scheme://")
            return True

    assert not Registry.get_format_class_for_file("scheme://something")
    registry.register(SchemeHandler)
    assert Registry.get_format_class_for_file("scheme://something") is SchemeHandler
    # Check dxtbx.load
    instance = dxtbx.load("scheme://something")
    assert isinstance(instance, SchemeHandler)


def test_multiple_scheme_handler(registry, tmp_path):
    _hits = []

    class SchemeHandler(Format):
        schemes = ["scheme", ""]

        @classmethod
        def understand(cls, endpoint):
            _hits.append(endpoint)
            return True

    registry.register(SchemeHandler)
    assert Registry.get_format_class_for_file("scheme://something") is SchemeHandler
    filename = tmp_path / "TEST_FILE"
    filename.touch()
    assert Registry.get_format_class_for_file(str(filename)) is SchemeHandler
    assert _hits == ["scheme://something", str(filename)]
