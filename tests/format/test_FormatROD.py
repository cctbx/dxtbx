"""Tests for dxbtx.format.FormatROD format classes."""

from __future__ import annotations

import pytest

from dxtbx.format.FormatROD import FormatROD, FormatROD_Arc
from dxtbx.model.experiment_list import ExperimentListFactory


def test_HyPix6000_TY6(dials_data):
    filename = (
        dials_data("image_examples") / "Hypix6000-monoclinic_lysozyme1_1_1.rodhypix.bz2"
    )

    assert FormatROD.understand(filename)
    expts = ExperimentListFactory.from_filenames([filename])
    assert len(expts) == 1

    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatROD

    gonio = imageset.get_goniometer()
    axes = gonio.get_names()
    assert list(axes) == ["PHI", "KAPPA=CHI", "OMEGA"]

    detector = imageset.get_detector()
    assert len(detector) == 1
    assert detector[0].get_gain() == 1.0
    assert detector[0].get_origin() == pytest.approx(
        (-44.6944, -36.9736, -23.2606), rel=1e-5
    )

    # TY6 decompression test
    img = imageset.get_raw_data(0)[0]
    assert img.focus() == (800, 775)
    assert img[123, 456] == 3


def test_EosS2_TY5(dials_data):
    filename = (
        dials_data("image_examples") / "EosS2-pre_SHR248e_CHCl3_pentane_1_1.img.bz2"
    )

    assert FormatROD.understand(filename)
    expts = ExperimentListFactory.from_filenames([filename])
    assert len(expts) == 1

    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatROD

    gonio = imageset.get_goniometer()
    axes = gonio.get_names()
    assert list(axes) == ["PHI", "KAPPA=CHI", "OMEGA"]

    detector = imageset.get_detector()
    assert len(detector) == 1
    # this CCD image has a non-unity gain
    assert detector[0].get_gain() == 100.0
    # if binning was not properly taken into account, the origin would break
    assert detector[0].get_origin() == pytest.approx(
        (-48.2605, -30.8897, -17.5460), rel=1e-5
    )

    # TY5 decompression test
    img = imageset.get_raw_data(0)[0]
    assert img.focus() == (512, 512)
    assert img[123, 456] == 174


def test_HyPixArc_TY6(dials_data):
    # HyPix-Arc is also tested in an end-to-end test `test_rigaku_hypix_arc_150()`
    # in dials/tests/algorithms/indexing/test_index.py.

    filename = dials_data("Rigaku_HyPix_Arc") / "CHNOS_1_2.rodhypix"

    assert FormatROD_Arc.understand(filename)
    expts = ExperimentListFactory.from_filenames([filename])
    assert len(expts) == 1

    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatROD_Arc

    # HyPix-Arc consists of three panels
    detector = imageset.get_detector()
    assert len(detector) == 3
