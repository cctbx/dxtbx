#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

import os
import socket

import pytest

collect_ignore = []


def dials_regression_path():
    """Return the absolute path to the dials_regression module as a string.
    This function is used directly by tests/test_regression_images.py"""
    try:
        import dials_regression as dr

        return os.path.abspath(os.path.dirname(dr.__file__))
    except ImportError:
        pass  # dials_regression not configured

    # Check if we are in a known location
    reference_copy = (
        "/dls/science/groups/scisoft/DIALS/repositories/git-reference/dials_regression"
    )
    if (
        os.name == "posix"
        and socket.gethostname().endswith(".diamond.ac.uk")
        and os.path.exists(reference_copy)
    ):
        return reference_copy


@pytest.fixture(scope="session")
def dials_regression():
    """Return the absolute path to the dials_regression module as a string.
    Skip the test if dials_regression is not available."""
    d_r = dials_regression_path()
    if d_r:
        return d_r
    pytest.skip("dials_regression required for this test")


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


@pytest.fixture
def run_in_tmpdir(tmpdir):
    """Shortcut to create a temporary directory and then run the test inside
    this directory."""
    tmpdir.chdir()
    return tmpdir
