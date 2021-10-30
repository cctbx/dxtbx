import copy
import pickle
import random

import pytest

from cctbx.eltbx import attenuation_coefficient
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.model import Beam, Detector, Panel, ParallaxCorrectedPxMmStrategy
from dxtbx.model.detector_helpers import (
    get_detector_projection_2d_axes,
    get_panel_projection_2d_from_axes,
    set_mosflm_beam_centre,
)
from dxtbx.model.experiment_list import ExperimentListFactory


def create_detector(offset):
    # Create the detector
    detector = Detector(
        Panel(
            "",  # Type
            "Panel",  # Name
            (10, 0, 0),  # Fast axis
            (0, 10, 0),  # Slow axis
            (0 + offset, 0 + offset, 200 - offset),
            # Origin
            (0.172, 0.172),  # Pixel size
            (512, 512),  # Image size
            (0, 1000),  # Trusted range
            0.1,  # Thickness
            "Si",  # Material
            identifier="123",
        )
    )
    return detector


def create_multipanel_detector(offset, ncols=3, nrows=2):
    reference_detector = create_detector(offset)
    reference_panel = reference_detector[0]
    panel_npix = reference_panel.get_image_size()
    panel_npix = int(panel_npix[0] / ncols), int(panel_npix[1] / nrows)

    multi_panel_detector = Detector()
    for i in range(ncols):
        for j in range(nrows):
            origin_pix = i * panel_npix[0], j * panel_npix[1]
            origin = reference_panel.get_pixel_lab_coord(origin_pix)
            new_panel = Panel(
                reference_panel.get_type(),
                reference_panel.get_name() + str(j + i * nrows),
                reference_panel.get_fast_axis(),
                reference_panel.get_slow_axis(),
                origin,
                reference_panel.get_pixel_size(),
                panel_npix,
                reference_panel.get_trusted_range(),
                reference_panel.get_thickness(),
                reference_panel.get_material(),
                identifier=reference_panel.get_identifier(),
            )
            multi_panel_detector.add_panel(new_panel)

    return multi_panel_detector


@pytest.fixture
def detector():
    return create_detector(offset=0)


def test_get_gain(detector):
    detector[0].set_gain(2.0)
    assert abs(detector[0].get_gain() - 2.0) < 1e-7


def test_get_pedestal(detector):
    detector[0].set_pedestal(2.0)
    assert abs(detector[0].get_pedestal() - 2.0) < 1e-7


def test_get_identifier(detector):
    detector[0].set_identifier("HELLO")
    assert detector[0].get_identifier() == "HELLO"
    detector2 = pickle.loads(pickle.dumps(detector))
    assert detector[0].get_identifier() == detector2[0].get_identifier()


def test_get_pixel_lab_coord(detector):
    eps = 1e-7

    # Check lab coordinates at the origin
    orig = detector[0].get_pixel_lab_coord((0, 0))
    dorig = abs(matrix.col(orig) - matrix.col(detector[0].get_origin()))
    assert dorig < eps

    # Check lab coordinate at the opposite corner
    corner = detector[0].get_pixel_lab_coord((512, 512))
    dcorner = abs(matrix.col(corner) - matrix.col(corner))
    assert dcorner < eps


def test_get_image_size_mm(detector):
    eps = 1e-7
    size = detector[0].get_image_size_mm()
    size2 = (512 * 0.172, 512 * 0.172)
    dsize = abs(matrix.col(size) - matrix.col(size2))
    assert dsize < eps


def test_is_value_in_trusted_range(detector):
    """Check values are either inside or outside trusted range."""
    assert detector[0].is_value_in_trusted_range(-1) is False
    assert detector[0].is_value_in_trusted_range(0) is True
    assert detector[0].is_value_in_trusted_range(999) is True
    assert detector[0].is_value_in_trusted_range(1000) is False


def test_is_coord_valid(detector):
    """Check points are either inside or outside detector range."""
    assert detector[0].is_coord_valid((-1, 256)) is False
    assert detector[0].is_coord_valid((256, 256)) is True
    assert detector[0].is_coord_valid((512, 256)) is False
    assert detector[0].is_coord_valid((256, -1)) is False
    assert detector[0].is_coord_valid((256, 256)) is True
    assert detector[0].is_coord_valid((256, 513)) is False


def test_pixel_to_millimeter_to_pixel(detector):
    eps = 1e-7

    # Pick some random pixels and check that px -> mm -> px give px == px
    w, h = detector[0].get_image_size()
    pixels = flex.vec2_double(
        (random.random() * w, random.random() * h) for i in range(100)
    )
    xy_mm = detector[0].pixel_to_millimeter(pixels)
    xy_px = detector[0].millimeter_to_pixel(xy_mm)
    assert approx_equal(xy_px, pixels, eps=eps)
    for xy in pixels:
        xy_mm = detector[0].pixel_to_millimeter(xy)
        xy_px = detector[0].millimeter_to_pixel(xy_mm)
        assert abs(matrix.col(xy_px) - matrix.col(xy)) < eps


def test_parallax_correction():
    # Attenuation length
    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(1) / 10.0
    t0 = 0.320

    # Create the detector
    detector = Detector(
        Panel(
            "",  # Type
            "",  # Name
            (10, 0, 0),  # Fast axis
            (0, 10, 0),  # Slow axis
            (0, 0, 200),  # Origin
            (0.172, 0.172),  # Pixel size
            (512, 512),  # Image size
            (0, 1000),  # Trusted range
            0.0,  # Thickness
            "",  # Material
            ParallaxCorrectedPxMmStrategy(mu, t0),
        )
    )
    for i in range(10000):
        mm = (random.uniform(-1000, 1000), random.uniform(-1000, 1000))
        px = detector[0].millimeter_to_pixel(mm)
        mm2 = detector[0].pixel_to_millimeter(px)
        assert abs(matrix.col(mm) - matrix.col(mm2)) < 1e-3


def test_get_names(detector):
    names = detector.get_names()
    assert len(names) == 1
    assert names[0] == "Panel"


def test_get_thickness(detector):
    for panel in detector:
        assert panel.get_thickness() == 0.1


def test_get_material(detector):
    for panel in detector:
        assert panel.get_material() == "Si"


def test_set_mosflm_beam_centre(detector):
    wavelength = 1.0
    panel = detector[0]
    detector_normal = matrix.col(panel.get_normal())
    _ = matrix.col(panel.get_origin())
    _ = matrix.col(panel.get_fast_axis())
    _ = matrix.col(panel.get_slow_axis())
    _ = panel.get_image_size_mm()

    s0 = (1.0 / wavelength) * detector_normal
    beam = Beam(-s0.normalize(), wavelength)

    beam_centre = matrix.col(panel.get_beam_centre(beam.get_s0()))
    origin_shift = matrix.col((1, 0.5))
    new_beam_centre = beam_centre + origin_shift

    new_mosflm_beam_centre = tuple(reversed(new_beam_centre))

    set_mosflm_beam_centre(detector, beam, new_mosflm_beam_centre)

    assert (
        matrix.col(panel.get_beam_centre(beam.get_s0()))
        - matrix.col(tuple(reversed(new_mosflm_beam_centre)))
    ).length() < 1e-6

    # test resolution methods
    beam = Beam(direction=(0, 0, 1), wavelength=1.0)
    d_min1 = detector.get_max_resolution(beam.get_s0())
    d_min2 = detector.get_max_inscribed_resolution(beam.get_s0())
    assert d_min1 < d_min2


def test_panel_mask():
    panel = Panel()
    panel.set_image_size((100, 100))
    panel.add_mask(40, 0, 60, 100)
    panel.add_mask(0, 40, 100, 60)
    panel.set_trusted_range((-1, 10))

    data = flex.double(flex.grid(100, 100))
    data[10, 10] = -1
    data[20, 20] = 10
    data[30, 30] = 100
    data[40, 40] = -10

    m1 = panel.get_untrusted_rectangle_mask()
    m2 = panel.get_trusted_range_mask(data)

    count = 0
    for j in range(100):
        for i in range(40, 60):
            assert m1[j, i] is False
            count += 1
    for i in range(100):
        for j in range(40, 60):
            if i >= 40 and i < 60:
                continue
            assert m1[j, i] is False
            count += 1
    assert m1.count(False) == count, "%d, %d" % (m1.count(False), count)

    assert m2.count(False) == 4
    assert m2[10, 10] is False
    assert m2[20, 20] is False
    assert m2[30, 30] is False
    assert m2[40, 40] is False


def test_equality():
    detector = create_detector(offset=0)

    # Create another detector with different origin
    # Equality operator on detector objects must find differences in origin
    detector_moved = create_detector(offset=100)
    assert detector != detector_moved

    # Equality operator on detector objects must identify identical detectors
    detector_moved_copy = create_detector(offset=100)
    assert detector_moved == detector_moved_copy


def test_panel_equality():
    panel = create_detector(offset=0)[0]
    panel2 = copy.deepcopy(panel)
    assert panel == panel2

    panel2.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(1, 1))
    assert panel != panel2

    panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(1, 1))
    assert panel == panel2


def test_get_detector_projection_2d_axes():
    # The function get_detector_projection_2d_axes should give the same results even if the
    # detector is rotated in the laboratory frame

    # Use a multipanel detector
    detector = create_multipanel_detector(offset=0)

    # Get 2D origin, fast and slow vectors for the detector
    o, f, s = get_detector_projection_2d_axes(detector)

    # Now rotate the detector by 30 degrees around an arbitrary axis
    h = detector.hierarchy()
    fast = matrix.col(h.get_fast_axis())
    slow = matrix.col(h.get_slow_axis())
    origin = matrix.col(h.get_origin())

    axis = matrix.col(flex.random_double_point_on_sphere())
    rot = axis.axis_and_angle_as_r3_rotation_matrix(30, deg=True)
    for panel in detector:
        fast = matrix.col(panel.get_fast_axis())
        slow = matrix.col(panel.get_slow_axis())
        origin = matrix.col(panel.get_origin())
        panel.set_frame(rot * fast, rot * slow, rot * origin)

    # Get 2D origin, fast and slow vectors for the rotated detector
    new_o, new_f, new_s = get_detector_projection_2d_axes(detector)

    for o1, o2 in zip(o, new_o):
        assert o1 == pytest.approx(o2)
    for f1, f2 in zip(f, new_f):
        assert f1 == pytest.approx(f2)
    for s1, s2 in zip(s, new_s):
        assert s1 == pytest.approx(s2)


def test_get_panel_projection_2d_from_axes(dials_data):

    # Get test data
    pytest.importorskip("h5py")
    filename = dials_data("image_examples", pathlib=True) / "dectris_eiger_master.h5"
    experiment = ExperimentListFactory.from_filenames([filename])[0]
    detector = experiment.detector

    # Get 2d axes
    origins, fast_axes, slow_axes = get_detector_projection_2d_axes(detector)

    # Get panel 0 data
    panel = detector[0]
    image_data = experiment.imageset.get_raw_data(0)[0]
    fast_axis = matrix.col(fast_axes[0] + (0,))
    slow_axis = matrix.col(slow_axes[0] + (0,))
    origin = matrix.col(origins[0] + (0,)) * 1e-3

    # Get 2d projection
    rotation, translation = get_panel_projection_2d_from_axes(
        panel, image_data, fast_axis, slow_axis, origin
    )

    expected_rotation = (
        1.0,
        0.0,
        0.0,
        1.0,
    )
    expected_translation = (1634.5, 1555.0)

    assert len(rotation) == len(expected_rotation)
    for rot, expected in zip(rotation, expected_rotation):
        assert rot == pytest.approx(expected)

    assert len(translation) == len(expected_translation)
    for t, expected in zip(translation, expected_translation):
        assert t == pytest.approx(expected)


def test_panel_get_projection_2d():

    detector = create_detector(offset=0)
    panel = detector[0]

    # Valid projection values
    valid_rotation = (1, 0, 0, 1)
    valid_translation = detector[0].get_image_size()

    # Test panel with no projections set explicitly has no 2d projection
    assert not panel.get_projection_2d()

    # Test panel with non-empty 2d projection has 2d projection
    panel.set_projection_2d(valid_rotation, valid_translation)
    assert panel.get_projection_2d()

    # Test values are set correctly
    rotation, translation = panel.get_projection_2d()
    assert rotation == valid_rotation
    assert translation == valid_translation


def test_detector_has_projection_2d():

    # Valid rotation value
    rotation = (1, 0, 0, 1)

    ## Test single panel detector
    detector = create_detector(offset=0)

    # Test detector with no projections set explicity has no 2d projection
    assert not detector.has_projection_2d()

    # Test detector with all panels having 2d projections gives a detector
    # with a 2d projection
    image_size = detector[0].get_image_size()
    for i in detector:
        i.set_projection_2d(rotation, (image_size[0], 0))

    assert detector.has_projection_2d()

    ## Test multipanel detector
    detector = create_multipanel_detector(offset=0)

    # Test detector with no projections set explicity has no 2d projection
    assert not detector.has_projection_2d()

    image_size = detector[0].get_image_size()
    for i in detector:
        # Test detector with only some panels with 2d projections gives a
        # detector without a 2d projection
        assert not detector.has_projection_2d()
        i.set_projection_2d(rotation, (image_size[0], 0))

    # Test detector with all panels having 2d projections gives a detector
    # with a 2d projection
    assert detector.has_projection_2d()
