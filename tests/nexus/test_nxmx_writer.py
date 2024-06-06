from __future__ import annotations

import glob

import h5py
import pytest

from libtbx.phil import parse
from libtbx.test_utils import approx_equal

from dxtbx.format.nexus import h5str
from dxtbx.format.nxmx_writer import NXmxWriter, phil_scope
from dxtbx.model.experiment_list import ExperimentListFactory

pytest.importorskip("dials")


def test_writer_jf16M(dials_data, tmpdir):
    h5path = (
        dials_data("lysozyme_JF16M_4img", pathlib=True)
        / "lyso009a_0087.JF07T32V01_master_4img.h5"
    )
    with h5py.File(h5path) as handle:
        output_file = tmpdir / "4img.h5"
        instrument_name = h5str(handle["entry/instrument/name"][()])
        source_name = h5str(handle["entry/source/name"][()])
        start_time = h5str(handle["entry/start_time"][()])
        end_time_estimated = h5str(handle["entry/end_time_estimated"][()])
        sample_name = h5str(handle["entry/sample/name"][()])
        params = phil_scope.fetch(
            parse(
                f"""
        output_file = {output_file}
        nexus_details {{
          instrument_name = {instrument_name}
          source_name = {source_name}
          start_time = {start_time}
          end_time_estimated = {end_time_estimated}
          sample_name = {sample_name}
        }}
        """
            )
        ).extract()

    expts1 = ExperimentListFactory.from_filenames([h5path])

    NXmxWriter(params)(expts1)

    expts2 = ExperimentListFactory.from_filenames([output_file])

    def recursive_test(pg1, pg2):
        assert pg1.is_similar_to(pg2)
        for c1, c2 in zip(pg1, pg2):
            recursive_test(c1, c2)

    for d1, d2 in zip(expts1.detectors(), expts2.detectors()):
        recursive_test(d1.hierarchy(), d2.hierarchy())

    for b1, b2 in zip(expts1.beams(), expts2.beams()):
        assert approx_equal(b1.get_wavelength(), b2.get_wavelength())
        assert approx_equal(b2.get_s0(), b2.get_s0())

    assert not (any(expts2.scans()))
    assert not (any(expts2.goniometers()))


def test_writer_x4wide(dials_data, tmpdir):
    cbfspath = dials_data("x4wide", pathlib=True) / "X4_wide_M1S4_2_000*.cbf"

    output_file = tmpdir / "x4wide_10img.h5"
    params = phil_scope.fetch(
        parse(
            f"""
    output_file = {output_file}
    nexus_details {{
      instrument_name = PILATUS
      source_name =  PILATUS 6M-F, S/N 60-0114-F
      start_time = 2013-02-08T13:23:40.833
      end_time_estimated = 2013-02-08T13:24:00.000
      sample_name = sample_x4wide
    }}
    """
        )
    ).extract()

    expts1 = ExperimentListFactory.from_filenames(glob.glob(str(cbfspath)))

    NXmxWriter(params)(expts1)

    expts2 = ExperimentListFactory.from_filenames([tmpdir / "x4wide_10img.h5"])

    def recursive_test(pg1, pg2):
        assert pg1.is_similar_to(pg2)
        for c1, c2 in zip(pg1, pg2):
            recursive_test(c1, c2)

    for d1, d2 in zip(expts1.detectors(), expts2.detectors()):
        recursive_test(d1.hierarchy(), d2.hierarchy())

    for b1, b2 in zip(expts1.beams(), expts2.beams()):
        assert approx_equal(b1.get_wavelength(), b2.get_wavelength())
        assert approx_equal(b2.get_s0(), b2.get_s0())

    for scan1, scan2 in zip(expts1.scans(), expts2.scans()):
        assert approx_equal(scan1.get_array_range(), scan2.get_array_range())
        assert approx_equal(scan1.get_oscillation(), scan2.get_oscillation())

    for gonio1, gonio2 in zip(expts1.goniometers(), expts2.goniometers()):
        assert approx_equal(gonio1.get_rotation_axis(), gonio2.get_rotation_axis())
