import math
import os
import pickle

import pytest

from scitbx import matrix
from scitbx.array_family import flex
from scitbx.math import principal_axes_of_inertia_2d

from dxtbx.masking import (
    GoniometerMaskerFactory,
    is_inside_polygon,
    mask_untrusted_polygon,
)
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model.goniometer import GoniometerFactory


def test_polygon():
    poly = flex.vec2_double(((0, 0), (1, 0), (1, 1), (0, 1)))
    assert is_inside_polygon(poly, 0, 0)
    assert is_inside_polygon(poly, 0.5, 0.5)
    assert not is_inside_polygon(poly, 1, 1.01)
    points = flex.vec2_double(((0.3, 0.8), (0.3, 1.5), (-8, 9), (0.00001, 0.9999)))
    assert list(is_inside_polygon(poly, points)) == [True, False, False, True]


@pytest.fixture
def kappa_goniometer():
    def _construct_goniometer(phi, kappa, omega):
        phi_axis = (1.0, 0.0, 0.0)
        kappa_axis = (0.914, 0.279, -0.297)
        omega_axis = (1.0, 0.0, 0.0)
        axes = flex.vec3_double((phi_axis, kappa_axis, omega_axis))
        angles = flex.double((phi, kappa, omega))
        names = flex.std_string(("GON_PHI", "GON_KAPPA", "GON_OMEGA"))
        return GoniometerFactory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    return _construct_goniometer


@pytest.fixture
def smargon_goniometer():
    def _construct_goniometer(phi, chi, omega):
        phi_axis = (1.0, 0.0, 0.0)
        chi_axis = (0, 0, -1)
        omega_axis = (1.0, 0.0, 0.0)
        axes = flex.vec3_double((phi_axis, chi_axis, omega_axis))
        angles = flex.double((phi, chi, omega))
        names = flex.std_string(("GON_PHI", "GON_CHI", "GON_OMEGA"))
        return GoniometerFactory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    return _construct_goniometer


@pytest.fixture
def pilatus_6M():
    def _construct_detector(distance):
        return DetectorFactory.simple(
            sensor="PAD",
            distance=distance,
            beam_centre=(216.87, 211.32),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(0.172, 0.172),
            image_size=(2463, 2527),
            trusted_range=(-1, 1e8),
        )

    return _construct_detector


@pytest.fixture(params=["cpp", "python"])
def kappa_goniometer_shadow_masker(request):
    def _construct_shadow_masker(goniometer):
        if request.param == "cpp":
            return GoniometerMaskerFactory.mini_kappa(goniometer)

        # Simple model of cone around goniometer phi axis
        # Exact values don't matter, only the ratio of height/radius
        height = 50  # mm
        radius = 20  # mm

        steps_per_degree = 1
        theta = (
            flex.double(range(360 * steps_per_degree))
            * math.pi
            / 180
            * 1
            / steps_per_degree
        )
        y = radius * flex.cos(-theta)
        z = radius * flex.sin(-theta)
        x = flex.double(theta.size(), height)

        coords = flex.vec3_double(zip(x, y, z))
        coords.insert(0, (0, 0, 0))

        return PyGoniometerShadowMasker(goniometer, coords, flex.size_t(len(coords), 0))

    return _construct_shadow_masker


@pytest.fixture(params=["cpp", "python"])
def smargon_shadow_masker(request):
    def _construct_shadow_masker(goniometer):
        if request.param == "python":
            return PySmarGonShadowMasker(goniometer)
        return GoniometerMaskerFactory.smargon(goniometer)

    return _construct_shadow_masker


@pytest.fixture
def dls_i23_experiment(dials_regression):
    experiments_file = os.path.join(
        dials_regression, "shadow_test_data/DLS_I23_Kappa/data_1_0400.cbf.gz"
    )
    experiments = ExperimentListFactory.from_filenames([experiments_file])
    return experiments[0]


@pytest.fixture
def dls_i23_kappa_shadow_masker():
    def _construct_shadow_masker(goniometer):
        return GoniometerMaskerFactory.dls_i23_kappa(goniometer)

    return _construct_shadow_masker


@pytest.fixture
def dls_i19_2_detector():
    def _construct_detector(distance):
        return DetectorFactory.simple(
            sensor="PAD",
            distance=distance,
            beam_centre=(41.20, 51.69),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(0.172, 0.172),
            image_size=(487, 619),
            trusted_range=(-1, 1e8),
        )

    return _construct_detector


@pytest.fixture
def dls_i19_2_goniometer():
    def _construct_goniometer(phi, kappa, omega):
        phi_axis = (1.0, 0.0, 0.0)
        kappa_axis = (0.642788, -0.766044, 0)
        omega_axis = (1.0, 0.0, 0.0)
        axes = flex.vec3_double((phi_axis, kappa_axis, omega_axis))
        angles = flex.double((phi, kappa, omega))
        names = flex.std_string(("GON_PHI", "GON_KAPPA", "GON_OMEGA"))
        return GoniometerFactory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    return _construct_goniometer


@pytest.fixture
def diamond_anvil_cell_masker():
    def _construct_shadow_masker(goniometer):
        return GoniometerMaskerFactory.diamond_anvil_cell(
            goniometer, cone_opening_angle=2 * 38 * math.pi / 180
        )

    return _construct_shadow_masker


def test_GoniometerShadowMasker_kappa_180_omega_0(
    kappa_goniometer, pilatus_6M, kappa_goniometer_shadow_masker
):
    goniometer = kappa_goniometer(phi=0, kappa=180, omega=0)
    detector = pilatus_6M(distance=170)
    masker = kappa_goniometer_shadow_masker(goniometer)

    scan_angle = 0  # omega = 0, shadow in top right corner of detector
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 361
    assert extrema[1] == pytest.approx(
        (43.60448791048147, 8.572923552543028, -30.416337975287743)
    )

    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow) == len(detector)
    assert len(shadow[0]) == 145

    mask = masker.get_mask(detector, scan_angle)
    assert len(mask) == len(detector)
    assert mask[0].all() == tuple(reversed(detector[0].get_image_size()))
    assert mask[0].count(True) == pytest.approx(5570865)

    obj = pickle.dumps(masker)
    masker2 = pickle.loads(obj)

    mask2 = masker2.get_mask(detector, scan_angle)
    assert len(mask2) == len(mask)
    assert mask2[0].all() == mask[0].all()
    assert mask2[0].count(True) == mask2[0].count(True)


def test_GoniometerShadowMasker_kappa_180_omega_m45(
    kappa_goniometer, pilatus_6M, kappa_goniometer_shadow_masker
):
    goniometer = kappa_goniometer(phi=0, kappa=180, omega=0)
    detector = pilatus_6M(distance=170)
    masker = kappa_goniometer_shadow_masker(goniometer)

    scan_angle = -45  # omega = -45, shadow on centre-right of detector
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 361
    assert extrema[1] == pytest.approx(
        (43.60448791048147, -15.445626462590827, -27.569571219784905)
    )
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow[0]) in (226, 227)
    # assert shadow[0][0] == pytest.approx((1623.5133257425446, 1254.0966982780835))
    mask = masker.get_mask(detector, scan_angle)
    assert mask[0].all() == tuple(reversed(detector[0].get_image_size()))
    assert mask[0].count(True) == pytest.approx(5467810)

    obj = pickle.dumps(masker)
    masker2 = pickle.loads(obj)

    mask2 = masker2.get_mask(detector, scan_angle)
    assert len(mask2) == len(mask)
    assert mask2[0].all() == mask[0].all()
    assert mask2[0].count(True) == mask2[0].count(True)


def test_GoniometerShadowMasker_kappa_180_omega_p45(
    kappa_goniometer, pilatus_6M, kappa_goniometer_shadow_masker
):
    goniometer = kappa_goniometer(phi=0, kappa=180, omega=0)
    detector = pilatus_6M(distance=170)
    masker = kappa_goniometer_shadow_masker(goniometer)

    scan_angle = 45  # omega = +45, no shadow on detector
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow[0]) == 0
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 361
    assert extrema[1] == pytest.approx(
        (43.60448791048147, 27.56957121978491, -15.445626462590822)
    )
    mask = masker.get_mask(detector, scan_angle)
    assert mask[0] is None or mask[0].count(False) == 0

    obj = pickle.dumps(masker)
    masker2 = pickle.loads(obj)

    mask2 = masker2.get_mask(detector, scan_angle)
    assert mask2[0] is None or mask2[0].count(False) == 0


def test_GoniometerShadowMasker_kappa_m70_omega_p100(
    kappa_goniometer, pilatus_6M, kappa_goniometer_shadow_masker
):
    # goniometer shadow does not intersect with detector panel
    goniometer = kappa_goniometer(phi=0, kappa=-70, omega=0)
    detector = pilatus_6M(distance=170)
    masker = kappa_goniometer_shadow_masker(goniometer)

    scan_angle = 100
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 361
    assert extrema[1] == pytest.approx(
        (42.318198019878935, 8.61724085750561, 32.17006801910824)
    )
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow[0]) == 0
    mask = masker.get_mask(detector, scan_angle)
    assert mask[0] is None or mask[0].count(False) == 0

    obj = pickle.dumps(masker)
    masker2 = pickle.loads(obj)

    mask2 = masker2.get_mask(detector, scan_angle)
    assert mask2[0] is None or mask2[0].count(False) == 0


def test_SmarGonShadowMasker_p48_c45_o95(
    smargon_goniometer, pilatus_6M, smargon_shadow_masker
):
    goniometer = smargon_goniometer(phi=48, chi=45, omega=100)
    detector = pilatus_6M(distance=170)
    masker = smargon_shadow_masker(goniometer)

    scan_angle = 100
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 82
    assert extrema[1] == pytest.approx(
        (22.106645739466337, 13.963679418895005, -22.47914302717216)
    )
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow[0]) in (15, 16)
    mask = masker.get_mask(detector, scan_angle)
    assert mask[0].count(True) == pytest.approx(5716721, 3e-5)

    obj = pickle.dumps(masker)
    masker2 = pickle.loads(obj)

    mask2 = masker2.get_mask(detector, scan_angle)
    assert len(mask2) == len(mask)
    assert mask2[0].all() == mask[0].all()
    assert mask2[0].count(True) == mask2[0].count(True)


def test_SmarGonShadowMasker_p0_c90_o50(
    smargon_goniometer, pilatus_6M, smargon_shadow_masker
):
    for phi in (0, 90, 180):
        # There should be no dependence of the shadow on phi
        goniometer = smargon_goniometer(phi=phi, chi=90, omega=50)
        detector = pilatus_6M(distance=214.75)
        masker = smargon_shadow_masker(goniometer)

        scan_angle = 50
        extrema = masker.extrema_at_scan_angle(scan_angle)
        assert len(extrema) == 82
        assert extrema[1] == pytest.approx(
            (-1.7364817766692957, -13.66792605230091, -31.609688838521162)
        )
        shadow = masker.project_extrema(detector, scan_angle)
        assert len(shadow[0]) == 11
        mask = masker.get_mask(detector, scan_angle)
        assert mask[0].count(True) == pytest.approx(4614865, 2e-4)

        obj = pickle.dumps(masker)
        masker2 = pickle.loads(obj)

        mask2 = masker2.get_mask(detector, scan_angle)
        assert len(mask2) == len(mask)
        assert mask2[0].all() == mask[0].all()
        assert mask2[0].count(True) == mask2[0].count(True)


def test_dls_i23_kappa(dls_i23_experiment, dls_i23_kappa_shadow_masker):
    for phi in (0, 90, 180):
        goniometer = dls_i23_experiment.goniometer
        goniometer.set_angles((0, -180, 0))
        detector = dls_i23_experiment.detector
        masker = dls_i23_kappa_shadow_masker(goniometer)

        scan_angle = 40
        extrema = masker.extrema_at_scan_angle(scan_angle)
        assert len(extrema) == 66
        assert extrema[1] == pytest.approx(
            (-18.68628039873836, -55.47017980234442, -98.76129898677576)
        )
        shadow = masker.project_extrema(detector, scan_angle)
        assert len(shadow) == len(detector)
        assert [len(s) for s in shadow] == [
            0,
            0,
            0,
            6,
            9,
            7,
            7,
            7,
            9,
            4,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        mask = masker.get_mask(detector, scan_angle)
        assert sum(m.count(False) for m in mask) == 1898843

        obj = pickle.dumps(masker)
        masker2 = pickle.loads(obj)

        mask2 = masker2.get_mask(detector, scan_angle)
        assert len(mask2) == len(mask)
        assert mask2[0].all() == mask[0].all()
        assert mask2[0].count(True) == mask2[0].count(True)


def test_dls_i19_2_diamond_anvil_cell(
    dls_i19_2_detector, dls_i19_2_goniometer, diamond_anvil_cell_masker
):
    detector = dls_i19_2_detector(distance=90.29)
    goniometer = dls_i19_2_goniometer(-35, 0, -90)
    masker = diamond_anvil_cell_masker(goniometer)
    scan_angle = -35

    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 721
    assert extrema[1] == pytest.approx(
        (7.812856265067172, 3.4202014332566875, -9.396926207859082)
    )
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow) == len(detector)
    assert len(shadow[0]) == 90
    mask = masker.get_mask(detector, scan_angle)
    assert sum(m.count(False) for m in mask) == 84352

    obj = pickle.dumps(masker)
    masker2 = pickle.loads(obj)

    mask2 = masker2.get_mask(detector, scan_angle)
    assert len(mask2) == len(mask)
    assert mask2[0].all() == mask[0].all()
    assert mask2[0].count(True) == mask2[0].count(True)


class PyGoniometerShadowMasker:
    def __init__(self, goniometer, extrema_at_datum, axis):
        # axis is an array of size_t the same size as extrema_at_datum,
        # where each element identifies the axis that that coordinate depends on
        self.goniometer = goniometer
        self._extrema_at_datum = extrema_at_datum
        self.axis = axis

    def extrema_at_scan_angle(self, scan_angle):
        axes = self.goniometer.get_axes()
        angles = self.goniometer.get_angles()
        scan_axis = self.goniometer.get_scan_axis()
        angles[scan_axis] = scan_angle
        extrema = self._extrema_at_datum.deep_copy()

        for i in range(len(axes)):
            sel = self.axis <= i
            rotation = matrix.col(axes[i]).axis_and_angle_as_r3_rotation_matrix(
                angles[i], deg=True
            )
            extrema.set_selected(sel, rotation.elems * extrema.select(sel))

        return extrema

    def project_extrema(self, detector, scan_angle):
        coords = self.extrema_at_scan_angle(scan_angle)
        shadow_boundary = []

        for p_id, p in enumerate(detector):
            # project coordinates onto panel plane
            a = p.get_D_matrix() * coords
            x, y, z = a.parts()
            valid = z > 0
            x.set_selected(valid, x.select(valid) / z.select(valid))
            y.set_selected(valid, y.select(valid) / z.select(valid))

            if valid.count(True) < 3:
                # no shadow projected onto this panel
                shadow_boundary.append(flex.vec2_double())
                continue

            # Compute convex hull of shadow points
            points = flex.vec2_double(x.select(valid), y.select(valid))
            shadow = flex.vec2_double(_convex_hull(points))
            shadow *= 1 / p.get_pixel_size()[0]

            shadow_orig = shadow.deep_copy()

            for i in (0, p.get_image_size()[0]):
                points = flex.vec2_double(
                    flex.double(p.get_image_size()[1], i),
                    flex.double_range(0, p.get_image_size()[1]),
                )
                inside = is_inside_polygon(shadow_orig, points)
                # only add those points needed to define vertices of shadow
                inside_isel = inside.iselection()
                outside_isel = (~inside).iselection()
                while inside_isel.size():
                    j = inside_isel[0]
                    shadow.append(points[j])
                    outside_isel = outside_isel.select(outside_isel > j)
                    if outside_isel.size() == 0:
                        shadow.append(points[inside_isel[-1]])
                        break
                    sel = inside_isel >= outside_isel[0]
                    if sel.count(True) == 0:
                        shadow.append(points[inside_isel[-1]])
                        break
                    inside_isel = inside_isel.select(sel)

            for i in (0, p.get_image_size()[1]):
                points = flex.vec2_double(
                    flex.double_range(0, p.get_image_size()[0]),
                    flex.double(p.get_image_size()[0], i),
                )
                inside = is_inside_polygon(shadow_orig, points)
                # only add those points needed to define vertices of shadow
                inside_isel = inside.iselection()
                outside_isel = (~inside).iselection()
                while inside_isel.size():
                    j = inside_isel[0]
                    shadow.append(points[j])
                    outside_isel = outside_isel.select(outside_isel > j)
                    if outside_isel.size() == 0:
                        shadow.append(points[inside_isel[-1]])
                        break
                    sel = inside_isel >= outside_isel[0]
                    if sel.count(True) == 0:
                        shadow.append(points[inside_isel[-1]])
                        break
                    inside_isel = inside_isel.select(sel)

            # Select only those vertices that are within the panel dimensions
            n_px = p.get_image_size()
            x, y = shadow.parts()
            valid = (x >= 0) & (x <= n_px[0]) & (y >= 0) & (y <= n_px[1])
            shadow = shadow.select(valid)

            # sort vertices clockwise from centre of mass
            com = principal_axes_of_inertia_2d(shadow).center_of_mass()
            sx, sy = shadow.parts()
            shadow = shadow.select(
                flex.sort_permutation(flex.atan2(sy - com[1], sx - com[0]))
            )

            shadow_boundary.append(shadow)

        return shadow_boundary

    def get_mask(self, detector, scan_angle):
        shadow_boundary = self.project_extrema(detector, scan_angle)

        mask = []
        for panel_id, panel in enumerate(detector):
            m = None
            if shadow_boundary[panel_id].size() > 3:
                m = flex.bool(flex.grid(reversed(panel.get_image_size())), True)
                mask_untrusted_polygon(m, shadow_boundary[panel_id])
            mask.append(m)
        return mask


# https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#Python
# https://github.com/thepracticaldev/orly-full-res/blob/master/copyingandpasting-big.png
def _convex_hull(points):
    """Computes the convex hull of a set of 2D points.

    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    """

    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list.
    return lower[:-1] + upper[:-1]


class PySmarGonShadowMasker(PyGoniometerShadowMasker):
    def __init__(self, goniometer):
        # FACE A: Sample holder
        #   Defined as semi-circle of radius r(A) = 10 mm (centred on PHI axis)
        #   with rectangle of size a(A) = 12.8 mm (x 20 mm)

        offsetA = 33.0
        # semi-circle for phi=-90 ... +90
        radiusA = 10.0
        phi = flex.double_range(-90, 100, step=10) * math.pi / 180
        x = flex.double(phi.size(), -offsetA)
        y = radiusA * flex.cos(phi)
        z = radiusA * flex.sin(phi)

        # corners of square
        sqdA = 12.8  # square depth
        nsteps = 10
        for i in range(nsteps + 1):
            for sign in (+1, -1):
                x.append(-offsetA)
                y.append(i * -sqdA / nsteps)
                z.append(sign * radiusA)
        x.append(-offsetA)
        y.append(-sqdA)
        z.append(0)

        self.faceA = flex.vec3_double(-x, -y, z)

        # FACE B: Lower arm
        sx = -28.50
        sy = -4.90
        sz = 8.50
        mx = -13.80
        my = -26.00
        nx = -27.50
        ny = -29.50
        px = -65.50
        py = -29.50
        self.faceB = flex.vec3_double(
            ((-sx, -sy, sz), (-mx, -my, 0), (-nx, -ny, 0), (-px, -py, 0))
        )

        # FACE E: Rim of sample holder
        #   Defined as circle of radius r(E) = 6 mm (centred on PHI axis) at an
        #   offset o(E) = 19 mm

        offsetE = 19.0
        radiusE = 6.0
        phi = flex.double_range(0, 360, step=15) * math.pi / 180
        x = flex.double(phi.size(), -offsetE)
        y = radiusE * flex.cos(phi)
        z = radiusE * flex.sin(phi)

        self.faceE = flex.vec3_double(-x, -y, z)

        extrema_at_datum = self.faceA.deep_copy()
        extrema_at_datum.extend(self.faceE)
        super().__init__(
            goniometer, extrema_at_datum, flex.size_t(extrema_at_datum.size(), 1)
        )

    def extrema_at_scan_angle(self, scan_angle):
        extrema = super().extrema_at_scan_angle(scan_angle)

        axes = self.goniometer.get_axes()
        angles = self.goniometer.get_angles()
        scan_axis = self.goniometer.get_scan_axis()
        angles[scan_axis] = scan_angle

        s = matrix.col(self.faceB[0])
        mx, my, _ = self.faceB[1]
        nx, ny, _ = self.faceB[2]
        px, py, _ = self.faceB[3]

        Rchi = matrix.col(axes[1]).axis_and_angle_as_r3_rotation_matrix(
            angles[1], deg=True
        )
        sk = Rchi * s
        sxk, syk, szk = sk.elems
        coords = flex.vec3_double(
            (
                (sxk, syk, 0),
                (sxk, syk, szk),
                (sxk + mx / 2, syk + my / 2, szk),
                (sxk + mx, syk + my, szk),
                (sxk + (mx + nx) / 2, syk + (my + ny) / 2, szk),
                (sxk + nx, syk + ny, szk),
                (sxk + (nx + px) / 2, syk + (ny + py) / 2, szk),
                (sxk + px, syk + py, szk),
                (sxk + px, syk + py, 0),
                (sxk + px, syk + py, -szk),
                (sxk + (nx + px) / 2, syk + (ny + py) / 2, -szk),
                (sxk + nx, syk + ny, -szk),
                (sxk + (mx + nx) / 2, syk + (my + ny) / 2, -szk),
                (sxk + mx, syk + my, -szk),
                (sxk + mx / 2, syk + my / 2, -szk),
                (sxk, syk, -szk),
            )
        )

        Romega = matrix.col(axes[2]).axis_and_angle_as_r3_rotation_matrix(
            angles[2], deg=True
        )
        coords = Romega.elems * coords
        extrema.extend(coords)

        return extrema
