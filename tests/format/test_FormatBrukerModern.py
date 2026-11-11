"""Tests for dxtbx.format.FormatBrukerModern"""

from __future__ import annotations

import bz2

import numpy as np
import pytest

from dxtbx.format.FormatBrukerModern import FormatBrukerModern
from dxtbx.model.experiment_list import ExperimentListFactory

HEADER_RECORD = 80


def _set_header_item(raw, key, value):
    """Replace the value of a header item in a raw Bruker image, keeping the
    80-character record length"""
    header_size = int(_get_header_item(raw, "HDRBLKS")) * 512
    tag = f"{key:<7}:".encode()
    for start in range(0, header_size, HEADER_RECORD):
        if raw[start : start + 8] == tag:
            record = (tag + value.encode()).ljust(HEADER_RECORD)
            assert len(record) == HEADER_RECORD
            return raw[:start] + record + raw[start + HEADER_RECORD :]
    raise KeyError(key)


def _get_header_item(raw, key):
    tag = f"{key:<7}:".encode()
    for start in range(0, len(raw), HEADER_RECORD):
        if raw[start : start + 8] == tag:
            return raw[start + 8 : start + HEADER_RECORD].decode().strip()
    raise KeyError(key)


def _pad16(b):
    return b + b"\0" * (-len(b) % 16)


def _mean(a):
    return float(np.mean(a.as_numpy_array()))


@pytest.fixture
def photon_ii_image(dials_data):
    return dials_data("image_examples") / "Bruker_PHOTON_II-dan_01_0001.sfrm.bz2"


@pytest.fixture
def photon_ii_raw(photon_ii_image):
    with bz2.open(photon_ii_image) as f:
        return f.read()


def _load(filename):
    assert FormatBrukerModern.understand(str(filename))
    expts = ExperimentListFactory.from_filenames([str(filename)])
    assert len(expts) == 1
    imageset = expts[0].imageset
    assert imageset.get_format_class() == FormatBrukerModern
    return imageset


def test_photon_ii_xray_panel_properties(photon_ii_image):
    """X-ray data should not have the electron diffraction-specific panel
    settings applied, and the baseline offset should be set as pedestal"""
    imageset = _load(photon_ii_image)
    assert imageset.get_beam().get_wavelength() > 0.1

    panel = imageset.get_detector()[0]
    assert panel.get_type() == "SENSOR_CCD"
    assert panel.get_mu() == 0.0
    assert panel.get_thickness() == 0.0
    assert panel.get_material() == ""
    assert panel.get_gain() == pytest.approx(137.6781 / 36.6)

    # NOVERFL[0] == -1 here, so the baseline NEXP[2] = 128 is stored in the
    # pixel values. This must be recorded as the pedestal.
    assert panel.get_pedestal() == 128.0
    raw = imageset.get_raw_data(0)[0]
    assert raw.focus() == (1024, 768)
    corrected = imageset.get_corrected_data(0)[0]
    assert _mean(corrected) == pytest.approx(
        (_mean(raw.as_double()) - 128.0) / panel.get_gain()
    )


def test_underflows_with_baseline_subtraction(photon_ii_raw, tmp_path):
    """Re-encode the image with baseline subtraction and an underflow table.
    Decoding should reproduce the original pixel values exactly."""
    original = _load_raw_data(photon_ii_raw, tmp_path / "orig_01_0001.sfrm")

    raw = photon_ii_raw
    header_size = int(_get_header_item(raw, "HDRBLKS")) * 512
    nrows = int(_get_header_item(raw, "NROWS").split()[0])
    ncols = int(_get_header_item(raw, "NCOLS").split()[0])
    noverfl = [int(e) for e in _get_header_item(raw, "NOVERFL").split()]
    assert _get_header_item(raw, "NPIXELB").split()[:2] == ["2", "1"]
    assert noverfl[:2] == [-1, 0]
    baseline = int(_get_header_item(raw, "NEXP").split()[2])
    assert baseline == 128

    npix = nrows * ncols
    data = np.frombuffer(raw, dtype="<u2", count=npix, offset=header_size)
    offset = header_size + npix * 2
    nbytes_4b = noverfl[2] * 4
    overflows_4b = np.frombuffer(raw, dtype="<u4", count=noverfl[2], offset=offset)
    tail = raw[offset + (-nbytes_4b % 16) + nbytes_4b :]

    # Subtract the baseline. Pixels flagged as 4-byte overflows keep the flag
    # value, with the baseline subtracted from the overflow table instead.
    is_overflow = data == 65535
    assert is_overflow.sum() == noverfl[2]
    full = data.astype(np.int64)
    subtracted = full - baseline
    is_underflow = (subtracted <= 0) & ~is_overflow
    num_underflows = int(is_underflow.sum())
    assert num_underflows > 0
    new_data = np.where(is_overflow, 65535, np.clip(subtracted, 0, None)).astype("<u2")
    # The underflow table stores the full, unsubtracted values
    underflow_table = full[is_underflow].astype(np.uint8)
    assert np.all(full[is_underflow] <= baseline)
    new_overflows_4b = (overflows_4b - baseline).astype("<u4")

    new_raw = _set_header_item(
        raw[:header_size], "NOVERFL", f"{num_underflows} 0 {noverfl[2]}"
    )
    new_raw += (
        new_data.tobytes()
        + _pad16(underflow_table.tobytes())
        + _pad16(new_overflows_4b.tobytes())
        + tail
    )

    reencoded = _load_raw_data(new_raw, tmp_path / "underflow_01_0001.sfrm")
    np.testing.assert_array_equal(reencoded, original)

    # The pedestal is the same whether or not the baseline was subtracted
    imageset = _load(tmp_path / "underflow_01_0001.sfrm")
    assert imageset.get_detector()[0].get_pedestal() == 128.0


def _load_raw_data(raw, path):
    path.write_bytes(raw)
    return _load(path).get_raw_data(0)[0].as_numpy_array()


def test_photon_iv_treated_as_photon(photon_ii_raw, photon_ii_image, tmp_path):
    """A Photon IV image should get the same models as the Photon II, not
    those intended for the EIGER-type ED-1 prototype"""
    dettype = _get_header_item(photon_ii_raw, "DETTYPE")
    assert dettype.startswith("CMOS-PHOTONII ")
    raw = _set_header_item(
        photon_ii_raw, "DETTYPE", dettype.replace("CMOS-PHOTONII", "CMOS-PHOTONIV")
    )
    path = tmp_path / "photon_iv_01_0001.sfrm"
    path.write_bytes(raw)

    photon_iv = _load(path)
    photon_ii = _load(photon_ii_image)

    assert photon_iv.get_detector()[0].get_type() == "SENSOR_CCD"
    assert photon_iv.get_detector()[0].get_pixel_size() == pytest.approx(
        (5.0 / 37.037037, 5.0 / 37.037037)
    )
    assert photon_iv.get_detector().is_similar_to(photon_ii.get_detector())
    assert photon_iv.get_goniometer().get_rotation_axis() == pytest.approx(
        photon_ii.get_goniometer().get_rotation_axis()
    )
    assert list(photon_iv.get_goniometer().get_axes().as_double()) == pytest.approx(
        list(photon_ii.get_goniometer().get_axes().as_double())
    )
