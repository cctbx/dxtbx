import functools
import math

from scitbx import matrix


def _all_equal(a, b):
    return all(x[0] == x[1] for x in zip(a, b))


def _all_approx_equal(a, b, tolerance):
    return all(abs(x[0] - x[1]) < tolerance for x in zip(a, b))


def sequence_diff(sequence1, sequence2, tolerance=None) -> str:
    b_diff = beam_diff
    d_diff = detector_diff
    g_diff = goniometer_diff
    s_diff = scan_diff
    if tolerance:
        b_diff = functools.partial(
            b_diff,
            wavelength_tolerance=tolerance.beam.wavelength,
            direction_tolerance=tolerance.beam.direction,
            polarization_normal_tolerance=tolerance.beam.polarization_normal,
            polarization_fraction_tolerance=tolerance.beam.polarization_fraction,
        )
        d_diff = functools.partial(
            d_diff,
            fast_axis_tolerance=tolerance.detector.fast_axis,
            slow_axis_tolerance=tolerance.detector.slow_axis,
            origin_tolerance=tolerance.detector.origin,
        )
        g_diff = functools.partial(
            g_diff,
            rotation_axis_tolerance=tolerance.goniometer.rotation_axis,
            fixed_rotation_tolerance=tolerance.goniometer.fixed_rotation,
            setting_rotation_tolerance=tolerance.goniometer.setting_rotation,
        )
        s_diff = functools.partial(s_diff, scan_tolerance=tolerance.scan.oscillation)
    output = (
        b_diff(sequence1.get_beam(), sequence2.get_beam()),
        d_diff(sequence1.get_detector(), sequence2.get_detector()),
        g_diff(sequence1.get_goniometer(), sequence2.get_goniometer()),
        s_diff(sequence1.get_scan(), sequence2.get_scan()),
    )

    return "\n".join(block for block in output if block)


def beam_diff(
    beam1,
    beam2,
    wavelength_tolerance=1e-6,
    direction_tolerance=1e-6,
    polarization_normal_tolerance=1e-6,
    polarization_fraction_tolerance=1e-6,
) -> str:
    aw = beam1.get_wavelength()
    bw = beam2.get_wavelength()
    ad = matrix.col(beam1.get_sample_to_source_direction())
    bd = matrix.col(beam2.get_sample_to_source_direction())
    an = matrix.col(beam1.get_polarization_normal())
    bn = matrix.col(beam2.get_polarization_normal())
    af = beam1.get_polarization_fraction()
    bf = beam2.get_polarization_fraction()
    text = []
    if abs(aw - bw) > wavelength_tolerance:
        text.append(f" Wavelength: {aw:f}, {bw:f}")
    if abs(ad.angle(bd)) > direction_tolerance:
        text.append(" Direction: {}, {}".format(tuple(ad), tuple(bd)))
    if abs(an.angle(bn)) > polarization_normal_tolerance:
        text.append(" Polarization Normal: {}, {}".format(tuple(an), tuple(bn)))
    if abs(af - bf) > polarization_fraction_tolerance:
        text.append(f" Polarization Fraction: {af}, {bf}")
    if len(text) > 0:
        text = ["Beam:"] + text
    return "\n".join(text)


def detector_diff(
    detector1,
    detector2,
    fast_axis_tolerance=1e-6,
    slow_axis_tolerance=1e-6,
    origin_tolerance=1e-6,
) -> str:
    text = []
    if len(detector1) != len(detector2):
        text.append("Num Panels: %d, %d" % (len(detector1), len(detector2)))
    for i, (aa, bb) in enumerate(zip(detector1, detector2)):
        a_image_size = aa.get_image_size()
        b_image_size = bb.get_image_size()
        a_pixel_size = aa.get_pixel_size()
        b_pixel_size = bb.get_pixel_size()
        a_trusted_range = aa.get_trusted_range()
        b_trusted_range = bb.get_trusted_range()
        a_fast = aa.get_fast_axis()
        b_fast = bb.get_fast_axis()
        a_slow = aa.get_slow_axis()
        b_slow = bb.get_slow_axis()
        a_origin = aa.get_origin()
        b_origin = bb.get_origin()
        temp_text = []
        if not _all_equal(a_image_size, b_image_size):
            temp_text.append(f"  Image size: {a_image_size}, {b_image_size}")
        if not _all_approx_equal(a_pixel_size, b_pixel_size, 1e-7):
            temp_text.append(f"  Pixel size: {a_pixel_size}, {b_pixel_size}")
        if not _all_approx_equal(a_trusted_range, b_trusted_range, 1e-7):
            temp_text.append(f"  Trusted Range: {a_trusted_range}, {b_trusted_range}")
        if not _all_approx_equal(a_fast, b_fast, fast_axis_tolerance):
            temp_text.append(f"  Fast axis: {a_fast}, {b_fast}")
        if not _all_approx_equal(a_slow, b_slow, slow_axis_tolerance):
            temp_text.append(f"  Slow axis: {a_slow}, {b_slow}")
        if not _all_approx_equal(a_origin, b_origin, origin_tolerance):
            temp_text.append(f"  Origin: {a_origin}, {b_origin}")
        if len(temp_text) > 0:
            text.append(" panel %d:" % i)
            text.extend(temp_text)
    if len(text) > 0:
        text = ["Detector:"] + text
    return "\n".join(text)


def goniometer_diff(
    goniometer1,
    goniometer2,
    rotation_axis_tolerance=1e-6,
    fixed_rotation_tolerance=1e-6,
    setting_rotation_tolerance=1e-6,
) -> str:
    a_axis = matrix.col(goniometer1.get_rotation_axis())
    b_axis = matrix.col(goniometer2.get_rotation_axis())
    a_fixed = goniometer1.get_fixed_rotation()
    b_fixed = goniometer2.get_fixed_rotation()
    a_setting = goniometer1.get_setting_rotation()
    b_setting = goniometer2.get_setting_rotation()
    text = []
    if abs(a_axis.angle(b_axis)) > rotation_axis_tolerance:
        text.append(" Rotation axis: {}, {}".format(tuple(a_axis), tuple(b_axis)))
    if not _all_approx_equal(a_fixed, b_fixed, fixed_rotation_tolerance):
        text.append(f" Fixed rotation: {a_fixed}, {b_fixed}")
    if not _all_approx_equal(a_setting, b_setting, setting_rotation_tolerance):
        text.append(f" Setting rotation: {a_setting}, {b_setting}")
    if len(text) > 0:
        text = ["Goniometer:"] + text
    return "\n".join(text)


def scan_diff(scan1, scan2, scan_tolerance=0.03) -> str:
    eps = scan_tolerance * abs(scan1.get_oscillation()[1])
    a_image_range = scan1.get_image_range()
    b_image_range = scan2.get_image_range()
    a_oscillation = scan1.get_oscillation()
    b_oscillation = scan2.get_oscillation()
    a_osc_range = scan1.get_oscillation_range()
    b_osc_range = scan2.get_oscillation_range()

    def mod_2pi(x):
        return x - 2 * math.pi * math.floor(x / (2 * math.pi))

    diff_2pi = abs(mod_2pi(a_osc_range[1]) - mod_2pi(b_osc_range[0]))
    diff_abs = abs(a_osc_range[1] - b_osc_range[0])
    text = []
    if not (a_image_range[1] + 1 == b_image_range[0]):
        text.append(f" Incompatible image range: {a_image_range}, {b_image_range}")
    if abs(a_oscillation[1] - b_oscillation[1]) > eps:
        text.append(f" Incompatible Oscillation: {a_oscillation}, {b_oscillation}")
    if min(diff_2pi, diff_abs) > eps * scan1.get_num_images():
        text.append(f" Incompatible Oscillation Range: {a_osc_range}, {b_osc_range}")
    if len(text) > 0:
        text = ["Scan:"] + text
    return "\n".join(text)
