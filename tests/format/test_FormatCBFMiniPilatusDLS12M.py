import calendar
import os

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

n_pixels_vertical_gaps = 195 * 7 * 4 * 24


@pytest.mark.parametrize(
    "timestamp,multi_panel,masked_pixel_count",
    (
        (None, True, 3053),
        (None, False, 3053 + n_pixels_vertical_gaps),
        (calendar.timegm((2019, 9, 3, 0, 0, 1)), True, 192957),
        (
            calendar.timegm((2019, 9, 3, 0, 0, 1)),
            False,
            192957 + n_pixels_vertical_gaps,
        ),
        (calendar.timegm((2019, 11, 26, 0, 0, 1)), True, 192948),
        (
            calendar.timegm((2019, 11, 26, 0, 0, 1)),
            False,
            192948 + n_pixels_vertical_gaps,
        ),
        (calendar.timegm((2020, 2, 21, 0, 0, 1)), True, 293699),
        (
            calendar.timegm((2020, 2, 21, 0, 0, 1)),
            False,
            293699 + n_pixels_vertical_gaps,
        ),
        (calendar.timegm((2020, 9, 8, 0, 0, 1)), True, 3053),
        (calendar.timegm((2020, 9, 8, 0, 0, 1)), False, 3053 + n_pixels_vertical_gaps),
    ),
)
def test_bad_pixel_mask(
    timestamp, multi_panel, masked_pixel_count, dials_regression, mocker
):
    if timestamp is not None:
        # Fool the format class into masking out bad modules
        mocked_timestamp = mocker.patch(
            "dxtbx.format.FormatCBFMiniPilatusDLS12M.get_pilatus_timestamp"
        )
        mocked_timestamp.return_value = timestamp

    image_path = os.path.join(
        dials_regression, "image_examples", "DLS_I23", "germ_13KeV_0001.cbf"
    )
    experiments = ExperimentListFactory.from_filenames(
        [image_path], format_kwargs={"multi_panel": multi_panel}
    )
    assert (
        sum(mask.count(False) for mask in experiments[0].imageset.get_mask(0))
        == masked_pixel_count
    )
