import os

import pytest

from dxtbx.model.scan_helpers import scan_helper_image_files
from dxtbx.model.sequence import SequenceFactory


@pytest.fixture
def image_test_dir(tmp_path):
    """Create a sequence of 20 'images' for scan tests to pick up"""
    for i in range(1, 21):
        image_file = tmp_path / f"image_{i:03d}.dat"
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
    """Test out the SequenceFactory."""

    template = "image_###.dat"

    xscans = [
        SequenceFactory.single_file(
            scan_helper_image_files.template_directory_index_to_image(
                template, image_test_dir, j + 1
            ),
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

    a = SequenceFactory.add(xscans[:10])
    b = SequenceFactory.add(xscans[10:])

    a + b

    filename = scan_helper_image_files.template_directory_index_to_image(
        template, image_test_dir, 1
    )

    assert len(SequenceFactory.search(filename)) == 20

    (a + b)[1:5]
    (a + b)[:10]
