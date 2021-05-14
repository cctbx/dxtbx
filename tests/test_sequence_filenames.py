import pytest

from dxtbx.sequence_filenames import template_regex


@pytest.mark.parametrize(
    "filename,template,digits",
    [
        ("foo_bar_001.img", "foo_bar_###.img", 1),
        ("foo_bar001.img", "foo_bar###.img", 1),
        ("foo_bar_1.8A_004.img", "foo_bar_1.8A_###.img", 4),
        ("foo_bar.001", "foo_bar.###", 1),
        ("foo_bar_002.img1000", "foo_bar_###.img1000", 2),
        ("foo_bar_00005.img", "foo_bar_#####.img", 5),
        ("image0010", "image####", 10),
    ],
)
def test_template_regex(filename, template, digits):
    assert template_regex(filename) == (template, digits)
