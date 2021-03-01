from dxtbx.datablock import BeamDiff, DetectorDiff, GoniometerDiff, ScanDiff


def sequence_diff(sequence1, sequence2, tolerance=None) -> str:
    if tolerance:
        b_diff = BeamDiff(
            wavelength_tolerance=tolerance.beam.wavelength,
            direction_tolerance=tolerance.beam.direction,
            polarization_normal_tolerance=tolerance.beam.polarization_normal,
            polarization_fraction_tolerance=tolerance.beam.polarization_fraction,
        )
        d_diff = DetectorDiff(
            fast_axis_tolerance=tolerance.detector.fast_axis,
            slow_axis_tolerance=tolerance.detector.slow_axis,
            origin_tolerance=tolerance.detector.origin,
        )
        g_diff = GoniometerDiff(
            rotation_axis_tolerance=tolerance.goniometer.rotation_axis,
            fixed_rotation_tolerance=tolerance.goniometer.fixed_rotation,
            setting_rotation_tolerance=tolerance.goniometer.setting_rotation,
        )
        s_diff = ScanDiff(scan_tolerance=tolerance.scan.oscillation)
    else:
        b_diff = BeamDiff()
        d_diff = DetectorDiff()
        g_diff = GoniometerDiff()
        s_diff = ScanDiff()

    return "\n".join(
        b_diff(sequence1.get_beam(), sequence2.get_beam())
        + d_diff(sequence1.get_detector(), sequence2.get_detector())
        + g_diff(sequence1.get_goniometer(), sequence2.get_goniometer())
        + s_diff(sequence1.get_scan(), sequence2.get_scan())
    )
