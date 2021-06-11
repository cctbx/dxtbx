"""
Helper methods for class for working with Pilatus images, for instance for
identifying the regions to be masked.
"""


import collections

_Detector = collections.namedtuple(
    "Detector",
    [
        "module_size_fast",
        "module_size_slow",
        "gap_fast",
        "gap_slow",
        "module_configs_fast",
    ],
)
setattr(
    _Detector,
    "all_widths",
    property(
        lambda self: [
            self.module_size_fast * j + self.gap_fast * (j - 1)
            for j in self.module_configs_fast
        ]
    ),
)

_DetectorDatabase = {
    "Pilatus": _Detector(
        module_size_fast=487,
        module_size_slow=195,
        gap_fast=7,
        gap_slow=17,
        module_configs_fast=(1, 3, 5),
    ),
    "Eiger X": _Detector(
        module_size_fast=1030,
        module_size_slow=514,
        gap_fast=10,
        gap_slow=37,
        module_configs_fast=(1, 2, 3, 4),
    ),
    "Eiger 2X": _Detector(
        module_size_fast=1028,
        module_size_slow=512,
        gap_fast=12,
        gap_slow=38,
        module_configs_fast=(1, 2, 3, 4),
    ),
}


def pilatus_6M_mask():
    """Hard coded mask regions for a Pilatus 6M instrument."""
    # FIX me, the parameters are listed here as f0, f1, s0, s1 but the prototype specifies f0, s0, f1, s1
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


def _get_pad_module_gap(xdetector, size=None):
    assert len(xdetector) == 1
    assert xdetector[0].get_type() == "SENSOR_PAD"

    if size is None:
        size = xdetector[0].get_image_size()

    for detector in _DetectorDatabase.values():
        if size[0] in detector.all_widths:
            return detector


def sensor_active_areas(xdetector):
    """Return the sensitive areas on the detector for pixel array detectors
    yes, does include hard coded magic numbers; returns [(x0, y0 x1, x1)]"""

    size = xdetector[0].get_image_size()
    detector = _get_pad_module_gap(xdetector)

    # iterate over these to produce a list of the TLC, BRC of the modules
    # as a list
    n_fast, remainder = divmod(size[0], detector.module_size_fast)
    assert (n_fast - 1) * detector.gap_fast == remainder

    n_slow, remainder = divmod(size[1], detector.module_size_slow)
    assert (n_slow - 1) * detector.gap_slow == remainder

    # Specify the dead areas between the modules, i.e. the rows and columns
    # where there are no active pixels

    panels = []

    for i_fast in range(n_fast):
        for i_slow in range(n_slow):
            panels.append(
                (
                    i_fast * detector.module_size_fast,
                    i_slow * detector.module_size_slow,
                    (i_fast + 1) * detector.module_size_fast,
                    (i_slow + 1) * detector.module_size_slow,
                )
            )

    return panels


def determine_pilatus_mask(xdetector):
    """Return an appropriate pixel mask for a Pilatus detector."""

    size = xdetector[0].get_image_size()

    # Hardcoded module size and gap size
    detector = _DetectorDatabase["Pilatus"]

    # Edge dead areas not included, only gaps between modules matter
    n_fast, remainder = divmod(size[0], detector.module_size_fast)
    assert (n_fast - 1) * detector.gap_fast == remainder

    n_slow, remainder = divmod(size[1], detector.module_size_slow)
    assert (n_slow - 1) * detector.gap_slow == remainder

    # Specify the dead areas between the modules, i.e. the rows and columns
    # where there are no active pixels
    mask = []
    for i_fast in range(n_fast - 1):
        mask.append(
            [
                (i_fast + 1) * detector.module_size_fast
                + i_fast * detector.gap_fast
                + 1,
                (i_fast + 1) * detector.module_size_fast
                + i_fast * detector.gap_fast
                + detector.gap_fast,
                1,
                size[1],
            ]
        )
    for i_slow in range(n_slow - 1):
        mask.append(
            [
                1,
                size[0],
                (i_slow + 1) * detector.module_size_slow
                + i_slow * detector.gap_slow
                + 1,
                (i_slow + 1) * detector.module_size_slow
                + i_slow * detector.gap_slow
                + detector.gap_slow,
            ]
        )

    return mask


def determine_eiger_mask(xdetector):
    """Return an appropriate pixel mask for an Eiger detector."""

    size = xdetector[0].get_image_size()

    detector = _get_pad_module_gap(xdetector)

    if not detector:
        size = tuple(reversed(xdetector[0].get_image_size()))
        detector = _get_pad_module_gap(xdetector, size=size)

    # Edge dead areas not included, only gaps between modules matter
    n_fast, remainder = divmod(size[0], detector.module_size_fast)
    assert (n_fast - 1) * detector.gap_fast == remainder

    n_slow, remainder = divmod(size[1], detector.module_size_slow)
    assert (n_slow - 1) * detector.gap_slow == remainder

    # Specify the dead areas between the modules, i.e. the rows and columns
    # where there are no active pixels
    mask = []
    for i_fast in range(n_fast - 1):
        mask.append(
            [
                (i_fast + 1) * detector.module_size_fast
                + i_fast * detector.gap_fast
                + 1,
                (i_fast + 1) * detector.module_size_fast
                + i_fast * detector.gap_fast
                + detector.gap_fast,
                1,
                size[1],
            ]
        )
    for i_slow in range(n_slow - 1):
        mask.append(
            [
                1,
                size[0],
                (i_slow + 1) * detector.module_size_slow
                + i_slow * detector.gap_slow
                + 1,
                (i_slow + 1) * detector.module_size_slow
                + i_slow * detector.gap_slow
                + detector.gap_slow,
            ]
        )

    return mask


def get_vendortype(xdetector):
    array = xdetector[0].get_image_size()
    if array == (2463, 2527):
        return "Pilatus-6M"
    elif array == (2463, 195):
        # special treatment of Pilatus-12M, treat as -6M for the viewer
        return "Pilatus-6M"
    elif array == (1475, 1679):
        return "Pilatus-2M"
    elif array == (487, 619):
        return "Pilatus-300K"
    return "Undetermined Pilatus size"


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
    elif array == (1030, 514):
        return "Eiger-500K"
    elif array == (4148, 4362):
        return "Eiger2-16M"
    elif array == (3108, 3262):
        return "Eiger2-9M"
    elif array == (2068, 2162):
        return "Eiger2-4M"
    elif array == (1028, 1062):
        return "Eiger2-1M"
    elif array == (1028, 512):
        return "Eiger2-500K"
    return "Undetermined EigerX size"
