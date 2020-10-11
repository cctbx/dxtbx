from __future__ import absolute_import, division, print_function

import os
from builtins import range

import pytest

from dxtbx.model.scan import ScanFactory
from dxtbx.model.scan_helpers import scan_helper_image_files


@pytest.fixture
def image_test_dir(tmp_path):
    """Create a sequence of 20 'images' for scan tests to pick up"""
    for i in range(1, 21):
        image_file = tmp_path / "image_{:03d}.dat".format(i)
        image_file.touch()
    return str(tmp_path)


def test_helper_image_files(image_test_dir):
    """Test the static methods in scan_helper_image_files."""
    assert scan_helper_image_files()

    template = "image_###.dat"

    assert (
        len(
            scan_helper_image_files.template_directory_to_indices(
                template, image_test_dir
            )
        )
        == 20
    )

    assert scan_helper_image_files.template_directory_index_to_image(
        template, image_test_dir, 1
    ) == os.path.join(image_test_dir, "image_001.dat")

    assert (
        scan_helper_image_files.template_index_to_image(template, 1) == "image_001.dat"
    )

    assert scan_helper_image_files.image_to_template_directory(
        os.path.join(image_test_dir, "image_001.dat")
    ) == (template, image_test_dir)

    assert scan_helper_image_files.image_to_index("image_001.dat") == 1

    assert scan_helper_image_files.image_to_template("image_001.dat") == "image_###.dat"

    assert scan_helper_image_files.image_to_index("image_6.8kev_1_001.cbf") == 1


def test_xScanFactory(image_test_dir):
    """Test out the ScanFactory."""

    template = "image_###.dat"

    xscans = [
        ScanFactory.single(
            scan_helper_image_files.template_directory_index_to_image(
                template, image_test_dir, j + 1
            ),
            None,
            1.0,
            18 + 0.5 * j,
            0.5,
            j,
        )
        for j in range(20)
    ]

    xscans.reverse()

    with pytest.raises(RuntimeError):
        print(sum(xscans[1:], xscans[0]))

    xscans.sort()

    sum(xscans[1:], xscans[0])

    a = ScanFactory.add(xscans[:10])
    b = ScanFactory.add(xscans[10:])

    a + b

    filename = scan_helper_image_files.template_directory_index_to_image(
        template, image_test_dir, 1
    )

    assert len(ScanFactory.search(filename)) == 20

    (a + b)[1:5]
    (a + b)[:10]
