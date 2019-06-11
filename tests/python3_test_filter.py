from __future__ import absolute_import, division, print_function

import pytest


class Python3TestFailureExpectationPlugin(object):
    def __init__(self, config):
        self.config = config
        self.known_failures_file = config.rootdir.join(".python3-test-failures")
        self.known_failures_output = config.rootdir.join(".python3-test-failures.new")
        if self.known_failures_file.check():
            self.known_failures = set(self.known_failures_file.readlines(cr=False))
        else:
            self.known_failures = set()
        self.revised_failures = set()

    def pytest_collectreport(self, report):
        if not report.nodeid:
            return
        passed = report.outcome in ("passed", "skipped")
        for item in report.result:
            if not passed:
                self.revised_failures.add(item.nodeid)

    def pytest_collection_modifyitems(self, session, config, items):
        for item in items:
            if item.nodeid in self.known_failures:
                item.add_marker(pytest.mark.xfail)

    def pytest_runtest_logreport(self, report):
        if report.failed:
            self.revised_failures.add(report.nodeid)
        if not report.passed and report.skipped and "xfail" in report.keywords:
            self.revised_failures.add(report.nodeid)

    def pytest_sessionfinish(self, session):
        self.known_failures_output.write("\n".join(sorted(self.revised_failures)))
