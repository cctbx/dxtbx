from __future__ import annotations

from dxtbx.model import Scan, ScanFactory


def test_scan_to_string_does_not_crash_on_empty_scan():
    print(Scan())


def test_scan_wrap_around_zero():
    filenames = [f"foobar_{n:03d}.cbf" for n in range(350, 370)]
    starts = [r % 360 for r in range(350, 370)]
    ends = [(r + 1) % 360 for r in range(350, 370)]
    scans = [
        ScanFactory.single_file(f, [0.1], s, e - s, 0)
        for f, s, e in zip(filenames, starts, ends)
    ]
    s0 = scans[0]
    for s in scans[1:]:
        assert s.get_oscillation()[1] == 1
        s0 += s
    assert s0.get_oscillation() == (350, 1)
    assert s0.get_image_range() == (350, 369)


def test_is_angle_valid():

    # Negative Scan

    scan = ScanFactory.make_scan(
        image_range=(1, 90),
        oscillation=(0, -1.0),
        exposure_times=0.1,
        epochs=range(90),
        deg=True,
    )

    expected_negative_range = [-i for i in range(91)]
    expected_positive_range = list(range(270, 360))
    total_range = list(range(-720, 720))

    for i in total_range:
        if i in expected_negative_range:
            assert scan.is_angle_valid(i)
        elif i in expected_positive_range or i % 360 in expected_positive_range:
            assert scan.is_angle_valid(i)
        else:
            assert scan.is_angle_valid(i) is False

    # Positive Scan

    scan = ScanFactory.make_scan(
        image_range=(1, 90),
        oscillation=(0, 1.0),
        exposure_times=0.1,
        epochs=range(90),
        deg=True,
    )
    expected_positive_range = list(range(91))
    expected_negative_range = list(range(-270, -360))
    for i in total_range:
        if i in expected_negative_range:
            assert scan.is_angle_valid(i)
        elif i in expected_positive_range or i % 360 in expected_positive_range:
            assert scan.is_angle_valid(i)
        else:
            assert scan.is_angle_valid(i) is False
