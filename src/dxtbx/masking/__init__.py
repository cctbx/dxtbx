import math

from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
from scitbx import matrix
from scitbx.array_family import flex

from dxtbx.model import MultiAxisGoniometer
from dxtbx_masking_ext import (
    GoniometerShadowMasker,
    SmarGonShadowMasker,
    is_inside_polygon,
    mask_untrusted_circle,
    mask_untrusted_polygon,
    mask_untrusted_rectangle,
)

__all__ = [
    "GoniometerShadowMasker",
    "SmarGonShadowMasker",
    "is_inside_polygon",
    "mask_untrusted_circle",
    "mask_untrusted_polygon",
    "mask_untrusted_rectangle",
]


class GoniometerMaskerFactory:
    @staticmethod
    def mini_kappa(goniometer, cone_opening_angle=43.60281897270362):
        """Construct a GoniometerShadowMasker for a mini-kappa goniometer.

        This is modelled a simple cone with the opening angle specified by
        `cone_opening_angle`.

        Args:
            goniometer (`dxtbx.model.Goniometer`): The goniometer instance.
            cone_opening_angle (float): The opening angle of the cone (in degrees).

        Returns:
            `dxtbx.masking.GoniometerShadowMasker`

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
            flex.double(range(360 * steps_per_degree))
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
            `dxtbx.masking.GoniometerShadowMasker`

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
            `dxtbx.masking.SmarGonShadowMasker`

        """
        return SmarGonShadowMasker(goniometer)

    @staticmethod
    def diamond_anvil_cell(goniometer, cone_opening_angle):
        radius_height_ratio = math.tan(1 / 2 * cone_opening_angle)
        height = 10  # mm
        radius = radius_height_ratio * height

        steps_per_degree = 1
        theta = (
            flex.double([list(range(360 * steps_per_degree))])
            * math.pi
            / 180
            * 1
            / steps_per_degree
        )
        x = radius * flex.cos(theta)  # x
        z = radius * flex.sin(theta)  # y
        y = flex.double(theta.size(), height)  # z

        coords = flex.vec3_double(zip(x, y, z))
        coords.extend(flex.vec3_double(zip(x, -y, z)))
        coords.insert(0, (0, 0, 0))

        return GoniometerShadowMasker(
            goniometer, coords, flex.size_t(len(coords), 0), True
        )
