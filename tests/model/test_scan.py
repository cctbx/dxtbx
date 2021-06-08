from __future__ import absolute_import, division, print_function

import dxtbx.model


def test_scan_to_string_does_not_crash_on_empty_scan():
    print(dxtbx.model.Scan())


def test_scan_wrap_around_zero():
    start = [r % 360 for r in range(350, 370)]
    end = [(r + 1) % 360 for r in range(350, 370)]
    scans = [
        dxtbx.model.scan.ScanFactory.single_file("n_123.x", [0.1], s, e - s, 0)
        for s, e in zip(start, end)
    ]
    for s in scans:
        assert s.get_oscillation()[1] == 1
