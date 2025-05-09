#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import annotations

import pytest

collect_ignore = []


def pytest_addoption(parser):
    """Add '--regression' options to pytest."""
    try:
        parser.addoption(
            "--regression",
            action="store_true",
            default=False,
            help="run (time-intensive) regression tests",
        )
    except ValueError:
        # Thrown in case the command line option is already defined
        pass


def pytest_collection_modifyitems(config, items):
    """Tests marked as regression are only run with --regression."""
    if not config.getoption("--regression"):
        skip_regression = pytest.mark.skip(reason="Test only runs with --regression")
        for item in items:
            if "regression" in item.keywords:
                item.add_marker(skip_regression)
