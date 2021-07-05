import itertools
import math
from operator import itemgetter

import numpy as np

from scitbx import matrix

try:
    import sklearn.cluster
except ImportError:
    sklearn = None


def read_xds_xparm(xds_xparm_file):
    """Parse the XDS XPARM file, which contains a description of the detector
    and experimental geometry, to a dictionary."""

    with open(xds_xparm_file) as fh:
        data = [float(i) for i in fh.read().split()]

    assert len(data) == 42

    starting_frame = int(data[0])
    phi_start, phi_width = data[1:3]
    axis = data[3:6]

    wavelength = data[6]
    beam = data[7:10]

    nx, ny = map(int, data[10:12])
    px, py = data[12:14]

    distance = data[14]
    ox, oy = data[15:17]

    x, y = data[17:20], data[20:23]
    normal = data[23:26]

    spacegroup = int(data[26])
    cell = data[27:33]

    a, b, c = data[33:36], data[36:39], data[39:42]

    results = {
        "starting_frame": starting_frame,
        "phi_start": phi_start,
        "phi_width": phi_width,
        "axis": axis,
        "wavelength": wavelength,
        "beam": beam,
        "nx": nx,
        "ny": ny,
        "px": px,
        "py": py,
        "distance": distance,
        "ox": ox,
        "oy": oy,
        "x": x,
        "y": y,
        "normal": normal,
        "spacegroup": spacegroup,
        "cell": cell,
        "a": a,
        "b": b,
        "c": c,
    }

    return results


def compute_frame_rotation(original, final):
    """Compute reference frame rotation to rotate from the original frame
    given by original = (x, y, z) to the to reference frame given by
    final = (_x, _y, _z). Returns M where M.x = _x etc."""

    x, y, z = original
    _x, _y, _z = final

    original_matrix = matrix.sqr(x.elems + y.elems + z.elems).transpose()
    assert (original_matrix.determinant() - 1.0) < 1.0e-7

    final_matrix = matrix.sqr(_x.elems + _y.elems + _z.elems).transpose()
    assert (final_matrix.determinant() - 1.0) < 1.0e-7

    # #1 rotate about x ^ (1, 0, 0) - if they are not coincident,
    # rotate about _x ^ _y if they are colinear but in opposite
    # directions

    if _x.angle(x) % math.pi:
        _ra_x = _x.cross(x)
        _a_x = _x.angle(x)
    elif math.fabs(_x.angle(x) - math.pi) < 1.0e-7:
        _ra_x = _x.cross(_y)
        _a_x = math.pi
    else:
        _ra_x = _x
        _a_x = 0.0

    _m_x = _ra_x.axis_and_angle_as_r3_rotation_matrix(-_a_x)

    # then rotate z to _z by rotating about _x (which is now coincident
    # with x)

    _ra_z = _x
    _a_z = _z.angle(_m_x * z)
    _m_z = _ra_z.axis_and_angle_as_r3_rotation_matrix(-_a_z)

    _m = _m_z * _m_x

    assert math.fabs(_m.determinant() - 1.0) < 1.0e-7

    return _m


def find_undefined_value(cbf_handle):
    """Given a cbf handle, get the value for the undefined pixel."""

    cbf_handle.find_category(b"array_intensities")
    cbf_handle.find_column(b"undefined_value")
    return cbf_handle.get_doublevalue()


def find_gain_value(cbf_handle):
    """Given a cbf handle, get the gain value."""
    try:
        cbf_handle.find_category(b"array_intensities")
        cbf_handle.find_column(b"gain")
    except Exception as e:
        if "CBF_NOTFOUND" not in str(e):
            raise
        return 1.0
    return cbf_handle.get_doublevalue()


class detector_helper_sensors:
    """A helper class which allows enumeration of detector sensor technologies
    which should help in identifying specific detectors when needed. These are
    currently limited to IMAGE_PLATE CCD PAD."""

    SENSOR_CCD = "SENSOR_CCD"
    SENSOR_PAD = "SENSOR_PAD"
    SENSOR_IMAGE_PLATE = "SENSOR_IMAGE_PLATE"
    SENSOR_UNKNOWN = "SENSOR_UNKNOWN"

    @staticmethod
    def check_sensor(sensor_type):
        if sensor_type in (
            detector_helper_sensors.SENSOR_CCD,
            detector_helper_sensors.SENSOR_PAD,
            detector_helper_sensors.SENSOR_IMAGE_PLATE,
            detector_helper_sensors.SENSOR_UNKNOWN,
        ):
            return True
        return False

    @staticmethod
    def all():
        return [
            detector_helper_sensors.SENSOR_CCD,
            detector_helper_sensors.SENSOR_PAD,
            detector_helper_sensors.SENSOR_IMAGE_PLATE,
        ]


def set_slow_fast_beam_centre_mm(detector, beam, beam_centre, panel_id=None):
    """detector and beam are dxtbx objects,
    beam_centre is a tuple of (slow, fast) mm coordinates.
    supports 2-theta offset detectors, assumes correct centre provided
    for 2-theta=0
    """
    beam_s, beam_f = beam_centre

    # Ensure panel_id is set
    us0 = matrix.col(beam.get_unit_s0())
    if panel_id is None:
        panel_id = detector.get_panel_intersection(us0)
        if panel_id < 0:
            panel_id = detector.get_panel_intersection(-us0)
            if panel_id < 0:
                panel_id = 0

    # Get data from the chosen panel
    panel = detector[panel_id]
    n = matrix.col(panel.get_normal())

    # Attempt to find the axis an angle of an applied 2theta shift
    cos_angle = n.cos_angle(us0)
    if cos_angle < 0:
        axi = us0.cross(-n)
        ang = us0.angle(-n)
    else:
        axi = us0.cross(n)
        ang = us0.angle(n)

    # Assume a 2theta offset if obliquity >= 5 deg
    two_theta = abs(ang) >= 5.0 * math.pi / 180.0

    # Undo 2theta shift
    if two_theta:
        R = axi.axis_and_angle_as_r3_rotation_matrix(ang)
        Rinv = R.inverse()
        try:
            h = detector.hierarchy()
            h.set_frame(
                fast_axis=Rinv * matrix.col(h.get_fast_axis()),
                slow_axis=Rinv * matrix.col(h.get_slow_axis()),
                origin=Rinv * matrix.col(h.get_origin()),
            )
        except AttributeError:
            for p in detector:
                p.set_frame(
                    fast_axis=Rinv * matrix.col(p.get_fast_axis()),
                    slow_axis=Rinv * matrix.col(p.get_slow_axis()),
                    origin=Rinv * matrix.col(p.get_origin()),
                )

    # Lab coord of desired beam centre
    if us0.accute_angle(n, deg=True) > 89.9:
        raise RuntimeError("Beam is in the plane of the detector panel")
    beam_dist = panel.get_directed_distance() / us0.dot(n)
    beam_centre_lab = beam_dist * us0

    # Lab coord of the current position where we want the beam centre
    intersection_lab = matrix.col(panel.get_lab_coord((beam_f, beam_s)))

    # If the detector has a hierarchy, just update the root note
    try:
        h = detector.hierarchy()
        translation = beam_centre_lab - intersection_lab
        new_origin = matrix.col(h.get_origin()) + translation
        h.set_frame(
            fast_axis=h.get_fast_axis(), slow_axis=h.get_slow_axis(), origin=new_origin
        )
    except AttributeError:
        # No hierarchy, update each panel instead by finding the offset of
        # its origin from the current position of the desired beam centre. Use
        # this to reposition the panel origin wrt the final beam centre
        for p in detector:
            origin = matrix.col(p.get_origin())
            offset = origin - intersection_lab
            new_origin = beam_centre_lab + offset
            p.set_frame(
                fast_axis=p.get_fast_axis(),
                slow_axis=p.get_slow_axis(),
                origin=new_origin,
            )

    # sanity check to make sure we have got the new beam centre correct
    new_beam_centre = detector[panel_id].get_bidirectional_ray_intersection(us0)
    assert (matrix.col(new_beam_centre) - matrix.col((beam_f, beam_s))).length() < 1e-4

    # Re-apply 2theta shift if required
    if two_theta:
        try:
            h = detector.hierarchy()
            h.set_frame(
                fast_axis=R * matrix.col(h.get_fast_axis()),
                slow_axis=R * matrix.col(h.get_slow_axis()),
                origin=R * matrix.col(h.get_origin()),
            )
        except AttributeError:
            for p in detector:
                p.set_frame(
                    fast_axis=R * matrix.col(p.get_fast_axis()),
                    slow_axis=R * matrix.col(p.get_slow_axis()),
                    origin=R * matrix.col(p.get_origin()),
                )

    return


def set_mosflm_beam_centre(detector, beam, mosflm_beam_centre):
    """detector and beam are dxtbx objects,
    mosflm_beam_centre is a tuple of mm coordinates.
    supports 2-theta offset detectors, assumes correct centre provided
    for 2-theta=0
    """
    return set_slow_fast_beam_centre_mm(detector, beam, mosflm_beam_centre)


def set_detector_distance(detector, distance):
    """
    Set detector origin from distance along normal
    """
    assert len(detector) == 1
    normal = matrix.col(detector[0].get_normal())
    origin = matrix.col(detector[0].get_origin())
    d = origin.dot(normal)
    x = origin - d * normal
    origin = distance * normal + x
    fast_axis = detector[0].get_fast_axis()
    slow_axis = detector[0].get_slow_axis()
    detector[0].set_frame(fast_axis, slow_axis, origin)


def project_2d(detector):
    """
    Project panel origin, fast and slow onto the best-fitting 2D plane.
    """

    # Extract panel vertices
    vertices = []
    for panel in detector:
        origin = matrix.col(panel.get_origin())
        fast = matrix.col(panel.get_fast_axis())
        slow = matrix.col(panel.get_slow_axis())
        panel_size = panel.get_image_size_mm()
        point1 = origin + panel_size[0] * fast
        point2 = origin + panel_size[1] * slow
        point3 = origin + panel_size[0] * fast + panel_size[1] * slow
        vertices.append(origin)
        vertices.append(point1)
        vertices.append(point2)
        vertices.append(point3)

    # Fit a plane by SVD. Modified from https://stackoverflow.com/a/18968498
    points = np.array(vertices).transpose()
    centre = points.mean(axis=1)
    r = points - centre[:, np.newaxis]  # centroid-to-vertex vectors
    inertia_tensor = np.matmul(r, r.T)
    u = np.linalg.svd(inertia_tensor)[0]
    normal = matrix.col(u[:, 2].tolist()).normalize()
    sum_vertices = matrix.col((0, 0, 0))
    for v in vertices:
        sum_vertices += v
    if normal.dot(sum_vertices) > 0:
        normal *= -1.0

    # For multi-panel detectors cluster fast, slow axes by DBSCAN to get a
    # consensus X, Y for the 2D plane
    clustered_axes = False
    if sklearn and len(detector) > 1:
        clustered_axes = True
        axes = []
        for panel in detector:
            axes.append(panel.get_fast_axis())
            axes.append(panel.get_slow_axis())
        clusters = sklearn.cluster.DBSCAN(eps=0.1, min_samples=2).fit_predict(axes)
        nclusters = max(clusters) + 1

        # Revert to single panel mode if clustering is unsucessful
        if nclusters < 2:
            clustered_axes = False

    if clustered_axes:
        summed_axes = [matrix.col((0, 0, 0)) for i in range(nclusters)]
        for axis, cluster in zip(axes, clusters):
            summed_axes[cluster] += matrix.col(axis)

        # Combine any two clusters approximately related by inversion
        all_checked = False
        while not all_checked:
            all_checked = True
            for i, j in itertools.combinations(range(len(summed_axes)), r=2):
                angle = summed_axes[i].angle(summed_axes[j], deg=True)
                if 175 < angle < 185:
                    summed_axes[i] = summed_axes[i] - summed_axes[j]
                    summed_axes = [e for k, e in enumerate(summed_axes) if k != j]
                    all_checked = False
                    break
        axes = summed_axes
    # For detectors with few panels or badly misaligned panels, align the
    # plane using the corners of the detector
    else:
        # Find the 4 points furthest from the centre (the detector corners)
        dists = np.linalg.norm(r, axis=0)
        threshold = sorted(dists)[-4]
        corners = np.where(dists >= threshold)[0].tolist()
        corners = [vertices[e] for e in corners]

        # Extract two spanning axes from these
        axes = [e - corners[0] for e in corners[1:]]
        lengths = [ax.length() for ax in axes]
        axes = [e for _, e in sorted(zip(lengths, axes), key=itemgetter(0))][0:2]

    # Choose whichever axis is closest to lab X for the plane's X
    labX = matrix.col((1, 0, 0))
    nearX = [ax.accute_angle(labX) for ax in axes]
    X = [ax for _, ax in sorted(zip(nearX, axes))][0]
    if X.dot(labX) < 0:
        X *= -1

    # Choose Y in the plane and X orthogonal
    Y = normal.cross(X).normalize()
    X = Y.cross(normal).normalize()

    # Project centre-shifted origins and fast, slow axes to the plane.
    origin_2d = []
    fast_2d = []
    slow_2d = []
    centre = matrix.col(centre)
    for panel in detector:
        origin = matrix.col(panel.get_origin())
        centre_to_origin = origin - centre
        fast = matrix.col(panel.get_fast_axis())
        slow = matrix.col(panel.get_slow_axis())

        origin_2d.append((centre_to_origin.dot(X), centre_to_origin.dot(Y)))
        fast_2d.append((fast.dot(X), fast.dot(Y)))
        slow_2d.append((slow.dot(X), slow.dot(Y)))

    return origin_2d, fast_2d, slow_2d
