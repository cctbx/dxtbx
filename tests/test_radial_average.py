from __future__ import annotations

from scitbx.array_family import flex

from dxtbx.ext import radial_average


def test_radial_average_keyword_arguments():
    data = flex.double(flex.grid(5, 5), 1)
    mask = flex.bool(flex.grid(5, 5), True)
    mask[0] = False
    sums = flex.double(10)
    sums_sq = flex.double(10)
    counts = flex.int(10)

    radial_average(
        data=data,
        mask=mask,
        beam_center=(2, 2),
        sums=sums,
        sums_sq=sums_sq,
        counts=counts,
        pixel_size=1.0,
        distance=100.0,
        upper_left=(0, 0),
        lower_right=(5, 5),
    )

    # One pixel is masked out
    assert flex.sum(counts) == 24
    assert flex.sum(sums) == 24
    assert flex.sum(sums_sq) == 24
