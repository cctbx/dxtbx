import dxtbx.model


def test_scan_to_string_does_not_crash_on_empty_scan():
    print(dxtbx.model.Scan())


def test_scan_wrap_around_zero():
    filenames = [f"foobar_{n:03d}.cbf" for n in range(350, 370)]
    starts = [r % 360 for r in range(350, 370)]
    ends = [(r + 1) % 360 for r in range(350, 370)]
    scans = [
        dxtbx.model.sequence.SequenceFactory.single_file(f, [0.1], s, e - s, 0)
        for f, s, e in zip(filenames, starts, ends)
    ]
    s0 = scans[0]
    for s in scans[1:]:
        assert s.get_oscillation()[1] == 1
        s0 += s
    assert s0.get_oscillation() == (350, 1)
    assert s0.get_image_range() == (350, 369)
