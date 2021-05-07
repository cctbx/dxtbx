from copy import deepcopy
from math import sqrt

import pytest
import six.moves.cPickle as pickle

from libtbx.phil import parse
from scitbx import matrix

from dxtbx.model import Beam, Detector
from dxtbx.model.detector import (
    DetectorFactory,
    ParallaxCorrectedPxMmStrategy,
    detector_phil_scope,
)


@pytest.fixture
def detector():
    detector = Detector()
    root = detector.hierarchy()
    root.set_name("D1")
    root.set_type("D")

    quad1 = root.add_group()
    quad1.set_name("Q1")
    quad1.set_type("Q")
    panel1 = quad1.add_panel()
    panel1.set_name("P1")
    panel1.set_type("P")
    panel2 = quad1.add_panel()
    panel2.set_name("P2")
    panel2.set_type("P")

    quad2 = root.add_group()
    quad2.set_name("Q2")
    quad2.set_type("Q")
    panel3 = quad2.add_panel()
    panel3.set_name("P3")
    panel3.set_type("P")
    panel4 = quad2.add_panel()
    panel4.set_name("P4")
    panel4.set_type("P")

    return detector


def test_flat(detector):
    """Test the flat hierarchy."""

    expected_types = ["P", "P", "P", "P"]
    expected_names = ["P1", "P2", "P3", "P4"]
    names = []
    types = []
    for p in detector:
        names.append(p.get_name())
        types.append(p.get_type())
    assert all(n == en for n, en in zip(names, expected_names))
    assert all(t == et for t, et in zip(types, expected_types))


def test_child_iteration(detector):
    # Iterate through the detector's children and check output
    expected_types = ["Q", "Q"]
    expected_names = ["Q1", "Q2"]
    names = []
    types = []
    for p in detector.hierarchy():
        names.append(p.get_name())
        types.append(p.get_type())
    assert all(n == en for n, en in zip(names, expected_names))
    assert all(t == et for t, et in zip(types, expected_types))


def test_child_iteration_in_reverse(detector):
    # Iterate through the detector's children in reverse and check output
    expected_types = ["Q", "Q"]
    expected_names = ["Q2", "Q1"]
    names = []
    types = []
    for p in reversed(detector.hierarchy()):
        names.append(p.get_name())
        types.append(p.get_type())
    assert all(n == en for n, en in zip(names, expected_names))
    assert all(t == et for t, et in zip(types, expected_types))


def test_child_indexing(detector):
    # Use an index to access the detector children
    assert len(detector.hierarchy()) == 2
    group = detector.hierarchy()[1]
    assert group.get_name() == "Q2" and group.get_type() == "Q"
    assert len(group) == 2
    panel = group[0]
    assert panel.get_name() == "P3" and panel.get_type() == "P"


def test_depth_first_iteration(detector):
    # Iterate through the tree pre-order and check output
    expected_types = ["D", "Q", "P", "P", "Q", "P", "P"]
    expected_names = ["D1", "Q1", "P1", "P2", "Q2", "P3", "P4"]
    names = []
    types = []
    for p in detector.iter_preorder():
        names.append(p.get_name())
        types.append(p.get_type())
    assert all(n == en for n, en in zip(names, expected_names))
    assert all(t == et for t, et in zip(types, expected_types))


def test_breadth_first_iteration(detector):
    # Iterate through the tree level-order and check output
    expected_types = ["D", "Q", "Q", "P", "P", "P", "P"]
    expected_names = ["D1", "Q1", "Q2", "P1", "P2", "P3", "P4"]
    names = []
    types = []
    for p in detector.iter_levelorder():
        names.append(p.get_name())
        types.append(p.get_type())
    assert all(n == en for n, en in zip(names, expected_names))
    assert all(t == et for t, et in zip(types, expected_types))


def test_panels_depth_first_iteration(detector):
    # Iterate through the panels in pre-order and check output
    expected_types = ["P", "P", "P", "P"]
    expected_names = ["P1", "P2", "P3", "P4"]
    names = []
    types = []
    for p in detector.iter_panels():
        names.append(p.get_name())
        types.append(p.get_type())
    assert all(n == en for n, en in zip(names, expected_names))
    assert all(t == et for t, et in zip(types, expected_types))


def test_get_valid_D_matrix(detector):
    """Setup the hierarchy of frames and check it's all consistent."""
    # Set a valid frame for the top level detector
    detector.hierarchy().set_local_frame(
        (1, 0, 0), (0, 1, 0), (0, 0, 100)  # Fast axis  # Slow axis
    )  # Origin

    # Check that all sub groups have the same frame and that we can get
    # a valid D matrix
    for obj in detector.iter_preorder():
        fast = matrix.col(obj.get_fast_axis())
        slow = matrix.col(obj.get_slow_axis())
        orig = matrix.col(obj.get_origin())
        assert abs(fast - matrix.col((1, 0, 0))) < 1e-7
        assert abs(slow - matrix.col((0, 1, 0))) < 1e-7
        assert abs(orig - matrix.col((0, 0, 100))) < 1e-7
        assert obj.get_D_matrix()

    # Get the quadrants and set their frames
    q1, q2 = detector.hierarchy().children()
    q1.set_local_frame(
        (1, 1, 0),  # Fast axis relative to detector frame
        (-1, 1, 0),  # Slow axis relative to detector frame
        (10, 10, 0),
    )  # Origin relative to detector frame
    q2.set_local_frame(
        (1, -1, 0),  # Fast axis relative to detector frame
        (1, 1, 0),  # Slow axis relative to detector frame
        (20, 20, 0),
    )  # Origin relative to detector frame

    # Get the panels and set their frames
    p1, p2 = q1.children()
    p1.set_local_frame(
        (1, -1, 0),  # Fast axis relative to q1 frame
        (1, 1, 0),  # Slow axis relative to q1 frame
        (5, 0, 10),
    )  # Origin relative to q1 frame
    p2.set_local_frame(
        (1, 1, 0),  # Fast axis relative to q1 frame
        (-1, 1, 0),  # Slow axis relative to q1 frame
        (0, 5, -10),
    )  # Origin relative to q1 frame

    # Get the panels and set their frames
    p3, p4 = q2.children()
    p3.set_local_frame(
        (1, -1, 0),  # Fast axis relative to q2 frame
        (1, 1, 0),  # Slow axis relative to q2 frame
        (0, 5, -10),
    )  # Origin relative to q2 frame
    p4.set_local_frame(
        (1, 1, 0),  # Fast axis relative to q2 frame
        (-1, 1, 0),  # Slow axis relative to q2 frame
        (5, 0, 10),
    )  # Origin relative to q2 frame

    # Test the panel coordinate systems

    eps = 1e-7
    p1_d0 = matrix.col((10.0 + sqrt(5.0 ** 2 / 2), 10.0 + sqrt(5.0 ** 2 / 2), 110))
    p2_d0 = matrix.col((10.0 - sqrt(5.0 ** 2 / 2), 10.0 + sqrt(5.0 ** 2 / 2), 90))
    p3_d0 = matrix.col((20.0 + sqrt(5.0 ** 2 / 2), 20.0 + sqrt(5.0 ** 2 / 2), 90))
    p4_d0 = matrix.col((20.0 + sqrt(5.0 ** 2 / 2), 20.0 - sqrt(5.0 ** 2 / 2), 110))
    p1_d1 = matrix.col((1, 0, 0))
    p2_d1 = matrix.col((0, 1, 0))
    p3_d1 = matrix.col((0, -1, 0))
    p4_d1 = matrix.col((1, 0, 0))
    p1_d2 = matrix.col((0, 1, 0))
    p2_d2 = matrix.col((-1, 0, 0))
    p3_d2 = matrix.col((1, 0, 0))
    p4_d2 = matrix.col((0, 1, 0))
    assert abs(matrix.col(p1.get_origin()) - p1_d0) < eps
    assert abs(matrix.col(p2.get_origin()) - p2_d0) < eps
    assert abs(matrix.col(p3.get_origin()) - p3_d0) < eps
    assert abs(matrix.col(p4.get_origin()) - p4_d0) < eps
    assert abs(matrix.col(p1.get_fast_axis()) - p1_d1) < eps
    assert abs(matrix.col(p2.get_fast_axis()) - p2_d1) < eps
    assert abs(matrix.col(p3.get_fast_axis()) - p3_d1) < eps
    assert abs(matrix.col(p4.get_fast_axis()) - p4_d1) < eps
    assert abs(matrix.col(p1.get_slow_axis()) - p1_d2) < eps
    assert abs(matrix.col(p2.get_slow_axis()) - p2_d2) < eps
    assert abs(matrix.col(p3.get_slow_axis()) - p3_d2) < eps
    assert abs(matrix.col(p4.get_slow_axis()) - p4_d2) < eps


def test_copy_and_reference(detector):
    # Get the detector hierarchy
    root = detector.hierarchy()

    # Get the panels in the hierarchy
    p1 = root[0][0]
    p2 = root[0][1]
    p3 = root[1][0]
    p4 = root[1][1]

    # Check panels are the same
    assert p1.is_(detector[0])
    assert p2.is_(detector[1])
    assert p3.is_(detector[2])
    assert p4.is_(detector[3])

    # Copy the detector
    new_detector = deepcopy(detector)

    # Check they're the same
    assert new_detector == detector

    # Add an offset to propagate
    root = new_detector.hierarchy()
    root.set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 10))

    # Get the panels in the hierarchy
    p1 = root[0][0]
    p2 = root[0][1]
    p3 = root[1][0]
    p4 = root[1][1]

    # Check panels are the same
    assert p1.is_(new_detector[0])
    assert p2.is_(new_detector[1])
    assert p3.is_(new_detector[2])
    assert p4.is_(new_detector[3])


def test_pickle(detector):
    # Get the detector hierarchy
    root = detector.hierarchy()

    # Get the panels in the hierarchy
    p1 = root[0][0]
    p2 = root[0][1]
    p3 = root[1][0]
    p4 = root[1][1]

    # Check panels are the same
    assert p1.is_(detector[0])
    assert p2.is_(detector[1])
    assert p3.is_(detector[2])
    assert p4.is_(detector[3])

    # Copy the detector
    new_detector = pickle.loads(pickle.dumps(detector))

    # Check they're the same
    assert new_detector == detector

    # Add an offset to propagate
    root = new_detector.hierarchy()
    root.set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 10))

    # Get the panels in the hierarchy
    p1 = root[0][0]
    p2 = root[0][1]
    p3 = root[1][0]
    p4 = root[1][1]

    # Check panels are the same
    assert p1.is_(new_detector[0])
    assert p2.is_(new_detector[1])
    assert p3.is_(new_detector[2])
    assert p4.is_(new_detector[3])


def test_from_phil():
    beam = Beam((0, 0, 1))

    params = detector_phil_scope.fetch(
        parse(
            """
    detector {
      panel {
        id = 0
        origin = (1, 1, 1)
        pixel_size = (0.001,0.001)
        image_size = (1000,1000)
        trusted_range = (-1, 1000)
        material = "Si"
        thickness = 0.01
        parallax_correction = True
      }
      panel {
        id = 1
        origin = (2, 2, 2)
        pixel_size = (0.001,0.001)
        image_size = (1000,1000)
        trusted_range = (-1, 1000)
      }
      panel {
        id = 2
        origin = (3, 3, 3)
        pixel_size = (0.001,0.001)
        image_size = (1000,1000)
        trusted_range = (-1, 1000)
      }
      panel {
        id = 3
        origin = (4, 4, 4)
        pixel_size = (0.001,0.001)
        image_size = (1000,1000)
        trusted_range = (-1, 1000)
      }

      hierarchy {
        name = "Root"
        origin = (100, 100, 100)
        group {
          id = 0
          origin = (10, 10, 10)
        }
        group {
          id = 0,0
          origin = (1, 1, 1)
          panel = 0
        }
        group {
          id = 0,1
          origin = (2, 2, 2)
          panel = 1
        }
        group {
          id = 1
          origin = (20, 20, 20)
        }
        group {
          id = 1,0
          origin = (1, 1, 1)
          panel = 2
        }
        group {
          id = 1,1
          origin = (2, 2, 2)
          panel = 3
        }
      }
    }
  """
        )
    ).extract()

    # Test create model
    d1 = DetectorFactory.from_phil(params, beam=beam)

    root = d1.hierarchy()

    # Check hierarchy origins
    assert root.get_origin() == (100, 100, 100)
    assert root[0].get_origin() == (110, 110, 110)
    assert root[1].get_origin() == (120, 120, 120)
    assert root[0][0].get_origin() == (111, 111, 111)
    assert root[0][1].get_origin() == (112, 112, 112)
    assert root[1][0].get_origin() == (121, 121, 121)
    assert root[1][1].get_origin() == (122, 122, 122)
    assert root[0][0][0].get_origin() == (112, 112, 112)
    assert root[0][1][0].get_origin() == (114, 114, 114)
    assert root[1][0][0].get_origin() == (124, 124, 124)
    assert root[1][1][0].get_origin() == (126, 126, 126)

    # Check panels are correct in hierarchy
    assert root[0][0][0].is_(d1[0])
    assert root[0][1][0].is_(d1[1])
    assert root[1][0][0].is_(d1[2])
    assert root[1][1][0].is_(d1[3])

    # Check panel attributes
    assert d1[0].get_image_size() == (1000, 1000)
    assert d1[0].get_pixel_size() == (0.001, 0.001)
    assert d1[0].get_trusted_range() == (-1, 1000)
    assert d1[0].get_material() == "Si"
    assert d1[0].get_thickness() == 0.01
    assert isinstance(d1[0].get_px_mm_strategy(), ParallaxCorrectedPxMmStrategy)

    assert d1[1].get_image_size() == (1000, 1000)
    assert d1[1].get_pixel_size() == (0.001, 0.001)
    assert d1[1].get_trusted_range() == (-1, 1000)

    assert d1[2].get_image_size() == (1000, 1000)
    assert d1[2].get_pixel_size() == (0.001, 0.001)
    assert d1[2].get_trusted_range() == (-1, 1000)

    assert d1[3].get_image_size() == (1000, 1000)
    assert d1[3].get_pixel_size() == (0.001, 0.001)
    assert d1[3].get_trusted_range() == (-1, 1000)

    params = detector_phil_scope.fetch(
        parse(
            """
    detector {
      panel {
        id = 0
        parallax_correction = False
      }
      panel {
        id = 1
        material = "Si"
        thickness = 0.01
        parallax_correction = True
      }

      hierarchy {
        name = "Root"
        origin = (200, 200, 200)
        group {
          id = 0
          origin = (20, 20, 20)
        }
        group {
          id = 0,0
          origin = (2, 2, 2)
        }
        group {
          id = 0,1
          origin = (3, 3, 3)
        }
        group {
          id = 1
          origin = (30, 30, 30)
        }
        group {
          id = 1,0
          origin = (2, 2, 2)
        }
        group {
          id = 1,1
          origin = (3, 3, 3)
        }
      }
    }
  """
        )
    ).extract()

    # Test overwrite model
    d2 = DetectorFactory.from_phil(params, reference=d1, beam=beam)

    root = d2.hierarchy()

    # Check hierarchy origins
    assert root.get_origin() == (200, 200, 200)
    assert root[0].get_origin() == (220, 220, 220)
    assert root[1].get_origin() == (230, 230, 230)
    assert root[0][0].get_origin() == (222, 222, 222)
    assert root[0][1].get_origin() == (223, 223, 223)
    assert root[1][0].get_origin() == (232, 232, 232)
    assert root[1][1].get_origin() == (233, 233, 233)
    assert root[0][0][0].get_origin() == (223, 223, 223)
    assert root[0][1][0].get_origin() == (225, 225, 225)
    assert root[1][0][0].get_origin() == (235, 235, 235)
    assert root[1][1][0].get_origin() == (237, 237, 237)

    # Check panels are correct in hierarchy
    assert root[0][0][0].is_(d2[0])
    assert root[0][1][0].is_(d2[1])
    assert root[1][0][0].is_(d2[2])
    assert root[1][1][0].is_(d2[3])

    # Check panel attributes
    assert not isinstance(d2[0].get_px_mm_strategy(), ParallaxCorrectedPxMmStrategy)
    assert d2[1].get_material() == "Si"
    assert d2[1].get_thickness() == 0.01
    assert isinstance(d2[1].get_px_mm_strategy(), ParallaxCorrectedPxMmStrategy)
