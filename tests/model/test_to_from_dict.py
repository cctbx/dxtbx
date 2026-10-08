from __future__ import annotations

import pytest

from scitbx.array_family import flex

from dxtbx.model import Beam, Detector, Goniometer, Scan


def test_beam():
    # Construct the object
    b1 = Beam(
        direction=(1, 2, 3),
        wavelength=1.1,
        divergence=0.01,
        sigma_divergence=0.01,
        polarization_normal=(4, 5, 6),
        polarization_fraction=0.9,
        flux=1.0,
        transmission=1.0,
    )

    # Create a dictionary and get the beam back
    d = b1.to_dict()
    b2 = Beam.from_dict(d)
    assert b2 == b1


def test_goniometer():
    # Construct the object
    g1 = Goniometer(
        rotation_axis=(1, 2, 3), fixed_rotation_matrix=(1, 2, 3, 4, 5, 6, 7, 8, 9)
    )

    # Create a dictionary and get the object back
    d = g1.to_dict()
    g2 = Goniometer.from_dict(d)
    assert g2 == g1


def test_scan():
    s1 = Scan(
        image_range=(1, 20),
        oscillation=(5.0, 0.1),
        exposure_times=flex.double(range(20)),
        epochs=flex.double(range(20, 40)),
    )
    s1.set_valid_image_ranges("0", [(1, 20)])

    d = s1.to_dict()
    s2 = Scan.from_dict(d)
    assert s1 == s2


def test_detector():
    d1 = Detector()
    p = d1.add_panel()
    p.set_name("p1")
    p.set_type("panel")
    p.set_pixel_size((0.1, 0.1))
    p.set_image_size((100, 100))
    p.set_trusted_range((0, 1000))
    p.set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))

    p = d1.add_panel()
    p.set_name("p2")
    p.set_type("panel")
    p.set_pixel_size((0.2, 0.2))
    p.set_image_size((200, 200))
    p.set_trusted_range((0, 2000))
    p.set_local_frame((0, 1, 0), (1, 0, 0), (0, 0, 1))

    root = d1.hierarchy()
    g = root.add_group()
    g.set_name("g1")
    g.set_type("group")
    g.set_local_frame((0, 1, 0), (1, 0, 0), (0, 0, 2))
    g.add_panel(d1[0])

    g = root.add_group()
    g.set_name("g2")
    g.set_type("group")
    g.set_local_frame((0, 1, 0), (1, 0, 0), (0, 0, 4))
    g.add_panel(d1[1])

    d = d1.to_dict()
    d2 = Detector.from_dict(d)
    assert len(d1) == len(d2)
    for p1, p2 in zip(d1, d2):
        assert p1 == p2
    assert d1.hierarchy() == d2.hierarchy()
    assert d1 == d2


def single_panel_detector():
    d = Detector()
    p = d.add_panel()
    p.set_name("p1")
    p.set_type("panel")
    p.set_pixel_size((0.1, 0.1))
    p.set_image_size((100, 100))
    p.set_trusted_range((0, 1000))
    p.set_frame((1, 0, 0), (0, -1, 0), (-5, 5, -100))
    return d


def test_single_panel_detector_is_not_written_with_a_hierarchy():
    """A hierarchy is meaningless for one panel. See dxtbx#472."""
    d1 = single_panel_detector()

    d = d1.to_dict()
    assert "hierarchy" not in d

    d2 = Detector.from_dict(d)
    assert d2 == d1
    assert not d2.has_hierarchy()


def test_detector_dict_written_before_dxtbx_472_still_loads():
    """Such a dict splits the geometry between the root node and the panel.

    The model itself still reads it back faithfully; it is
    DetectorFactory.from_dict that folds the two together.
    """
    d = single_panel_detector().to_dict()
    d["hierarchy"] = {
        "fast_axis": (1.0, 0.0, 0.0),
        "slow_axis": (0.0, 1.0, 0.0),
        "origin": (0.0, 0.0, -100.0),
        "children": [{"panel": 0}],
    }

    detector = Detector.from_dict(d)
    assert len(detector) == 1
    assert detector.has_hierarchy()
    # The panel dict holds the local frame, so the root is added on top of it
    assert detector[0].get_origin() == pytest.approx((-5, 5, -200))
