from __future__ import annotations

import shutil

import pytest

from dxtbx.sequence_filenames import template_image_range, template_regex


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
        ("foo_123_1_1.rodhypix", "foo_123_1_#.rodhypix", 1),  # Rigaku-style
    ],
)
def test_template_regex(filename, template, digits):
    assert template_regex(filename) == (template, digits)


def test_template_image_range(dials_data):
    template = str(dials_data("insulin", pathlib=True) / "insulin_1_###.img")
    assert template_image_range(template) == (1, 45)


def test_template_image_range_non_zero_padded(dials_data, tmp_path):
    images = sorted(dials_data("insulin", pathlib=True).glob("insulin_1_0[0-1]*"))
    # symlink if possible, copy if necessary
    for i, image in enumerate(images):
        try:
            (tmp_path / f"insulin_1_{i + 1}.img").symlink_to(image)
        except OSError:
            shutil.copy(image, (tmp_path / f"insulin_1_{i + 1}.img"))

    template = str(tmp_path / "insulin_1_#.img")
    assert template_image_range(template) == (1, 19)
