from __future__ import absolute_import, division

#!/usr/bin/env python
# FormatPilatusHelpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Helper methods for class for working with Pilatus images, for instance for
# identifying the regions to be masked.


def pilatus_6M_mask():
    """Hard coded mask regions for a Pilatus 6M instrument."""
    # FIX me, the paramters are listed here as f0, f1, s0, s1 but the prototype specifies f0, s0, f1, s1
    return [
        [488, 494, 1, 2527],
        [982, 988, 1, 2527],
        [1476, 1482, 1, 2527],
        [1970, 1976, 1, 2527],
        [1, 2463, 196, 212],
        [1, 2463, 408, 424],
        [1, 2463, 620, 636],
        [1, 2463, 832, 848],
        [1, 2463, 1044, 1060],
        [1, 2463, 1256, 1272],
        [1, 2463, 1468, 1484],
        [1, 2463, 1680, 1696],
        [1, 2463, 1892, 1908],
        [1, 2463, 2104, 2120],
        [1, 2463, 2316, 2332],
    ]


def pilatus_2M_mask():
    """Hard coded mask regions for a Pilatus 2M detector."""

    return [
        [488, 494, 1, 1679],
        [982, 988, 1, 1679],
        [1, 1475, 196, 212],
        [1, 1475, 408, 424],
        [1, 1475, 620, 636],
        [1, 1475, 832, 848],
        [1, 1475, 1044, 1060],
        [1, 1475, 1256, 1272],
        [1, 1475, 1468, 1484],
    ]


def pilatus_300K_mask():
    """Hard coded mask regions for a Pilatus 300K instrument."""

    return [[1, 487, 196, 212], [1, 487, 408, 424]]


def get_pad_module_gap(xdetector, size=None):
    assert xdetector[0].get_type() == "SENSOR_PAD"
    assert len(xdetector) == 1

    if size is None:
        size = xdetector[0].get_image_size()

    pilatus_width = 487
    pilatus_height = 195
    pilatus_gap = 7

    pilatus_sizes = [pilatus_width * j + pilatus_gap * (j - 1) for j in (1, 3, 5)]

    eiger_x_width = 1030
    eiger_x_height = 514
    eiger_x_gap = 10

    eiger_2x_width = 1028
    eiger_2x_height = 512
    eiger_2x_gap = 12

    eiger_x_sizes = [eiger_x_width * j + eiger_x_gap * (j - 1) for j in (1, 2, 3, 4)]
    eiger_2x_sizes = [eiger_2x_width * j + eiger_2x_gap * (j - 1) for j in (1, 2, 3, 4)]

    if size[0] in pilatus_sizes:
        module_size_fast, module_size_slow = (487, 195)
        gap_size_fast, gap_size_slow = (7, 17)
    elif size[0] in eiger_x_sizes:
        module_size_fast, module_size_slow = (1030, 514)
        gap_size_fast, gap_size_slow = (10, 37)
    elif size[0] in eiger_2x_sizes:
        module_size_fast, module_size_slow = (1028, 512)
        gap_size_fast, gap_size_slow = (12, 38)
    else:
        return None

    return (module_size_fast, module_size_slow, gap_size_fast, gap_size_slow)


def sensor_active_areas(xdetector):
    """Return the sensitive areas on the detector for pixel array detectors
    yes, does include hard coded magic numbers; returns [(x0, y0 x1, x1)]"""

    _ = get_pad_module_gap(xdetector)
    module_size_fast, module_size_slow, gap_size_fast, gap_size_slow = _

    # iterate over these to produce a list of the TLC, BRC of the modules
    # as a list
    n_fast, remainder = divmod(size[0], module_size_fast)
    assert (n_fast - 1) * gap_size_fast == remainder

    n_slow, remainder = divmod(size[1], module_size_slow)
    assert (n_slow - 1) * gap_size_slow == remainder

    # Specify the dead areas between the modules, i.e. the rows and columns
    # where there are no active pixels

    panels = []

    for i_fast in range(n_fast):
        for i_slow in range(n_slow):
            panels.append(
                (
                    i_fast * module_size_fast,
                    i_slow * module_size_slow,
                    (i_fast + 1) * module_size_fast,
                    (i_slow + 1) * module_size_slow,
                )
            )

    return panels


def determine_pilatus_mask(xdetector):
    """Return an appropriate pixel mask for a Pilatus detector."""

    size = xdetector[0].get_image_size()

    # Hardcoded module size and gap size
    module_size_fast, module_size_slow = (487, 195)
    gap_size_fast, gap_size_slow = (7, 17)

    # Edge dead areas not included, only gaps between modules matter
    n_fast, remainder = divmod(size[0], module_size_fast)
    assert (n_fast - 1) * gap_size_fast == remainder

    n_slow, remainder = divmod(size[1], module_size_slow)
    assert (n_slow - 1) * gap_size_slow == remainder

    # Specify the dead areas between the modules, i.e. the rows and columns
    # where there are no active pixels
    mask = []
    for i_fast in range(n_fast - 1):
        mask.append(
            [
                (i_fast + 1) * module_size_fast + i_fast * gap_size_fast + 1,
                (i_fast + 1) * module_size_fast
                + i_fast * gap_size_fast
                + gap_size_fast,
                1,
                size[1],
            ]
        )
    for i_slow in range(n_slow - 1):
        mask.append(
            [
                1,
                size[0],
                (i_slow + 1) * module_size_slow + i_slow * gap_size_slow + 1,
                (i_slow + 1) * module_size_slow
                + i_slow * gap_size_slow
                + gap_size_slow,
            ]
        )

    return mask


def determine_eiger_mask(xdetector):
    """Return an appropriate pixel mask for an Eiger detector."""

    size = xdetector[0].get_image_size()

    _ = get_pad_module_gap(xdetector)

    if _ is None:
        size = tuple(reversed(xdetector[0].get_image_size()))
        _ = get_pad_module_gap(xdetector, size=size)

    module_size_fast, module_size_slow, gap_size_fast, gap_size_slow = _

    # Edge dead areas not included, only gaps between modules matter
    n_fast, remainder = divmod(size[0], module_size_fast)
    assert (n_fast - 1) * gap_size_fast == remainder

    n_slow, remainder = divmod(size[1], module_size_slow)
    assert (n_slow - 1) * gap_size_slow == remainder

    # Specify the dead areas between the modules, i.e. the rows and columns
    # where there are no active pixels
    mask = []
    for i_fast in range(n_fast - 1):
        mask.append(
            [
                (i_fast + 1) * module_size_fast + i_fast * gap_size_fast + 1,
                (i_fast + 1) * module_size_fast
                + i_fast * gap_size_fast
                + gap_size_fast,
                1,
                size[1],
            ]
        )
    for i_slow in range(n_slow - 1):
        mask.append(
            [
                1,
                size[0],
                (i_slow + 1) * module_size_slow + i_slow * gap_size_slow + 1,
                (i_slow + 1) * module_size_slow
                + i_slow * gap_size_slow
                + gap_size_slow,
            ]
        )

    return mask


def get_vendortype(xdetector):
    array = xdetector[0].get_image_size()
    if array == (2463, 2527):
        return "Pilatus-6M"
    elif array == (2463, 195):
        return (
            "Pilatus-6M"
        )  # special treatment of Pilatus-12M, treat as -6M for the viewer
    elif array == (1475, 1679):
        return "Pilatus-2M"
    elif array == (487, 619):
        return "Pilatus-300K"
    return "Undetermined Pilatus size"


# FIXME this contains hard coded information and does not apply to 2X
def get_vendortype_eiger(xdetector):
    array = xdetector[0].get_image_size()
    # print array,
    if array == (4150, 4371):
        return "Eiger-16M"
    elif array == (3110, 3269):
        return "Eiger-9M"
    elif array == (2070, 2167):
        return "Eiger-4M"
    elif array == (1030, 1065):
        return "Eiger-1M"
    return "Undetermined EigerX size"
