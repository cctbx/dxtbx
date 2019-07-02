from __future__ import absolute_import, division, print_function

from dxtbx_util_masking_ext import (  # noqa: F401, exported symbol
    mask_untrusted_rectangle,
    mask_untrusted_circle,
    mask_untrusted_polygon,
    is_inside_polygon,
)

import math

from scitbx import matrix
from scitbx.math import principal_axes_of_inertia_2d
from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame

from dxtbx.model import MultiAxisGoniometer

from scitbx.array_family import flex
from dxtbx_util_masking_ext import GoniometerShadowMasker  # noqa: F401
from dxtbx_util_masking_ext import SmarGonShadowMasker  # noqa: F401


class PyGoniometerShadowMasker(object):
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
        for panel_id in range(len(detector)):
            m = None
            if shadow_boundary[panel_id].size() > 3:
                m = flex.bool(
                    flex.grid(reversed(detector[panel_id].get_image_size())), True
                )
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


class GoniometerMaskerFactory(object):
    @staticmethod
    def mini_kappa(goniometer, cone_opening_angle=43.60281897270362):
        """Construct a GoniometerShadowMasker for a mini-kappa goniometer.

        This is modelled a simple cone with the opening angle specified by
        `cone_opening_angle`.

        Args:
            goniometer (`dxtbx.model.Goniometer`): The goniometer instance.
            cone_opening_angle (float): The opening angle of the cone (in degrees).

        Returns:
            `dxtbx.util.masking.GoniometerShadowMasker`

        """
        assert isinstance(goniometer, MultiAxisGoniometer)
        assert len(goniometer.get_axes()) == 3
        # Simple model of cone around goniometer phi axis
        # Exact values don't matter, only the ratio of height/radius
        height = 50  # mm

        radius_height_ratio = math.tan(1 / 2 * cone_opening_angle * math.pi / 180)
        radius = radius_height_ratio * height

        steps_per_degree = 1
        theta = (
            flex.double([range(360 * steps_per_degree)])
            * math.pi
            / 180
            * 1
            / steps_per_degree
        )
        y = radius * flex.cos(-theta)  # x
        z = radius * flex.sin(-theta)  # y
        x = flex.double(theta.size(), height)  # z

        coords = flex.vec3_double(zip(x, y, z))
        coords.insert(0, (0, 0, 0))

        return GoniometerShadowMasker(goniometer, coords, flex.size_t(len(coords), 0))

    @staticmethod
    def dls_i23_kappa(goniometer):
        """Construct a GoniometerShadowMasker for the DLS I23 Kappa goniometer.

        Args:
            goniometer (`dxtbx.model.Goniometer`): The goniometer instance.

        Returns:
            `dxtbx.util.masking.GoniometerShadowMasker`

        """
        coords = flex.vec3_double(((0, 0, 0),))

        alpha = flex.double_range(0, 190, step=10) * math.pi / 180
        r = flex.double(alpha.size(), 40)
        x = flex.double(r.size(), 107.61)
        y = -r * flex.sin(alpha)
        z = -r * flex.cos(alpha)
        coords.extend(flex.vec3_double(x, y, z))

        coords.extend(
            flex.vec3_double(
                (
                    # fixed
                    (107.49, 7.84, 39.49),
                    (107.39, 15.69, 38.97),
                    (107.27, 23.53, 38.46),
                    (107.16, 31.37, 37.94),
                    (101.76, 33.99, 36.25),
                    (96.37, 36.63, 34.56),
                    (90.98, 39.25, 33.00),
                    (85.58, 41.88, 31.18),
                    (80.89, 47.06, 31.00),
                    (76.55, 51.51, 31.03),
                    (72.90, 55.04, 31.18),
                    (66.86, 60.46, 31.67),
                    (62.10, 64.41, 32.25),
                )
            )
        )

        alpha = flex.double_range(180, 370, step=10) * math.pi / 180
        r = flex.double(alpha.size(), 33)
        x = flex.sqrt(flex.pow2(r * flex.sin(alpha)) + 89.02 ** 2) * flex.cos(
            (50 * math.pi / 180) - flex.atan(r / 89.02 * flex.sin(alpha))
        )
        y = flex.sqrt(flex.pow2(r * flex.sin(alpha)) + 89.02 ** 2) * flex.sin(
            (50 * math.pi / 180) - flex.atan(r / 89.02 * flex.sin(alpha))
        )
        z = -r * flex.cos(alpha)
        coords.extend(flex.vec3_double(x, y, z))

        coords.extend(
            flex.vec3_double(
                (
                    # fixed
                    (62.10, 64.41, -32.25),
                    (66.86, 60.46, -31.67),
                    (72.90, 55.04, -31.18),
                    (76.55, 51.51, -31.03),
                    (80.89, 47.06, -31.00),
                    (85.58, 41.88, -31.18),
                    (90.98, 39.25, -33.00),
                    (96.37, 36.63, -34.56),
                    (101.76, 33.99, -36.25),
                    (107.16, 31.37, -37.94),
                    (107.27, 23.53, -38.46),
                    (107.39, 15.69, -38.97),
                    (107.49, 7.84, -39.49),
                    (107.61, 0.00, -40.00),
                )
            )
        )

        # I23 end station coordinate system:
        #   X-axis: positive direction is facing away from the storage ring (from
        #           sample towards goniometer)
        #   Y-axis: positive direction is vertically up
        #   Z-axis: positive direction is in the direction of the beam (from
        #           sample towards detector)
        #   K-axis (kappa): at an angle of +50 degrees from the X-axis
        #   K & phi rotation axes: clockwise rotation is positive (right hand
        #           thumb rule)
        #   Omega-axis: along the X-axis; clockwise rotation is positive

        # End station x-axis is parallel to ImgCIF x-axis
        # End station z-axis points in opposite direction to ImgCIF definition
        # (ImgCIF: The Z-axis is derived from the source axis which goes from
        # the sample to the source)
        # Consequently end station y-axis (to complete set following right hand
        # rule) points in opposite direction to ImgCIF y-axis.
        # Kappa arm aligned with -y in ImgCIF convention

        R = align_reference_frame(
            matrix.col((1, 0, 0)),
            matrix.col((1, 0, 0)),
            matrix.col((0, 1, 0)),
            matrix.col((0, -1, 0)),
        )
        coords = R.elems * coords

        return GoniometerShadowMasker(goniometer, coords, flex.size_t(len(coords), 1))

    @staticmethod
    def smargon(goniometer):
        """Construct a SmarGonShadowMasker for the SmarGon goniometer.

        Args:
            goniometer (`dxtbx.model.Goniometer`): The goniometer instance.

        Returns:
            `dxtbx.util.masking.SmarGonShadowMasker`

        """
        return SmarGonShadowMasker(goniometer)
