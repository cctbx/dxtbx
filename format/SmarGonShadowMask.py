from __future__ import absolute_import, division

from builtins import range
import math

from dials.util.masking import GoniometerShadowMaskGenerator


class SmarGonShadowMaskGenerator(GoniometerShadowMaskGenerator):
    def __init__(self, goniometer):
        from scitbx.array_family import flex

        self.goniometer = goniometer

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

        self.faceA = flex.vec3_double(x, y, z)

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
            ((sx, sy, sz), (mx, my, 0), (nx, ny, 0), (px, py, 0))
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

        self.faceE = flex.vec3_double(x, y, z)

    def extrema_at_scan_angle(self, scan_angle):
        from scitbx.array_family import flex

        # Align end station coordinate system with ImgCIF coordinate system
        from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
        from scitbx import matrix

        R = align_reference_frame(
            matrix.col((-1, 0, 0)),
            matrix.col((1, 0, 0)),
            matrix.col((0, -1, 0)),
            matrix.col((0, 1, 0)),
        )
        faceA = R.elems * self.faceA
        faceE = R.elems * self.faceE

        axes = self.goniometer.get_axes()
        angles = self.goniometer.get_angles()
        scan_axis = self.goniometer.get_scan_axis()
        angles[scan_axis] = scan_angle

        extrema = flex.vec3_double()

        for coords in (faceA, faceE):
            coords = coords.deep_copy()
            for i, axis in enumerate(axes):
                if i == 0:
                    continue  # shadow doesn't change with phi setting
                sel = flex.bool(len(coords), True)
                rotation = matrix.col(axis).axis_and_angle_as_r3_rotation_matrix(
                    angles[i], deg=True
                )
                coords.set_selected(sel, rotation.elems * coords.select(sel))
            extrema.extend(coords)

        s = matrix.col(self.faceB[0])
        mx, my, _ = self.faceB[1]
        nx, ny, _ = self.faceB[2]
        px, py, _ = self.faceB[3]

        Rchi = (R.inverse() * matrix.col(axes[1])).axis_and_angle_as_r3_rotation_matrix(
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

        coords = R.elems * coords
        Romega = matrix.col(axes[2]).axis_and_angle_as_r3_rotation_matrix(
            angles[2], deg=True
        )
        coords = Romega.elems * coords
        extrema.extend(coords)

        return extrema
