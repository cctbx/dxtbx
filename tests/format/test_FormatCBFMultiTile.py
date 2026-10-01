"""Round-trip a multi-panel detector whose panels have different pixel counts
through the full-CBF writer and the FormatCBFMultiTile(Hierarchy) readers."""

from __future__ import annotations

import numpy as np
import pytest

from scitbx.array_family import flex

import dxtbx
from dxtbx.format.cbf_writer import FullCBFWriter
from dxtbx.format.FormatCBFMultiTile import FormatCBFMultiTileStill
from dxtbx.format.FormatCBFMultiTileHierarchy import FormatCBFMultiTileHierarchyStill
from dxtbx.imageset import ImageSet, ImageSetData, MemReader
from dxtbx.model import Beam, Detector, Panel


class _InMemoryPanels:
    def __init__(self, panels):
        self.panels = tuple(panels)

    def get_raw_data(self):
        return self.panels

    def get_mask(self, goniometer=None):
        return tuple(flex.bool(flex.grid(p.focus()), True) for p in self.panels)


class _NamedMemReader(MemReader):
    def paths(self):
        return ["in_memory_%d" % i for i in range(len(self._images))]


def _make_imageset(sizes_fast_slow):
    beam = Beam()
    beam.set_unit_s0((0, 0, -1))
    beam.set_wavelength(1.0)
    det = Detector()
    data = []
    off = 0.0
    for i, (nf, ns) in enumerate(sizes_fast_slow):
        p = Panel()
        p.set_name("p%d" % i)
        p.set_image_size((nf, ns))
        p.set_pixel_size((0.1, 0.1))
        p.set_frame((1, 0, 0), (0, -1, 0), (-5 + off, 5, -100))
        p.set_trusted_range((-1, 1e6))
        det.add_panel(p)
        off += nf * 0.1 + 1.0
        arr = flex.int(flex.grid(ns, nf), 0)
        arr[0, 0] = 1000 + i  # tag each panel
        data.append(arr)
    reader = _NamedMemReader([_InMemoryPanels(data)])
    reader.format_class = _InMemoryPanels
    imageset = ImageSet(ImageSetData(reader, None))
    imageset.set_beam(beam)
    imageset.set_detector(det)
    return imageset, det, data


def test_multitile_cbf_panels_with_different_sizes(tmp_path):
    sizes = [(40, 30), (20, 30), (40, 10)]
    imageset, det, data = _make_imageset(sizes)

    filename = str(tmp_path / "irregular.cbf")
    writer = FullCBFWriter(imageset=imageset)
    cbf = writer.get_cbf_handle(index=0, header_only=True)
    writer.add_data_to_cbf(cbf, data=tuple(data))
    writer.write_cbf(filename, cbf=cbf)

    # default registry choice is the hierarchy reader: check models and data
    fmt = dxtbx.load(filename)
    assert isinstance(fmt, FormatCBFMultiTileHierarchyStill)
    detector = fmt.get_detector()
    raw = fmt.get_raw_data()
    assert len(detector) == len(raw) == len(sizes)
    for panel, expected, arr, truth in zip(detector, sizes, raw, det):
        assert panel.get_image_size() == expected
        assert arr.focus() == (expected[1], expected[0])
        assert np.allclose(panel.get_origin(), truth.get_origin())
    assert [arr[0, 0] for arr in raw] == [1000 + i for i in range(len(sizes))]

    # the plain (non-hierarchy) multi-tile reader must also get per-panel sizes
    plain = FormatCBFMultiTileStill(filename)
    assert [p.get_image_size() for p in plain.get_detector()] == sizes


def test_multitile_cbf_rotation_round_trip(tmp_path):
    """Goniometer and scan written to a multi-panel full CBF come back as a
    rotation sequence with the same models."""
    from dxtbx.model import Goniometer, Scan
    from dxtbx.model.experiment_list import ExperimentListFactory

    sizes = [(40, 30), (20, 30)]
    imageset, det, data = _make_imageset(sizes)
    gonio = Goniometer((0, 1, 0))
    n = 3
    scan = Scan(
        (1, n),
        (10.0, 0.25),
        exposure_times=flex.double(n, 0.5),
        epochs=flex.double([1.7e9 + i for i in range(n)]),
    )
    filenames = []
    for i in range(n):
        filename = str(tmp_path / ("rot_%04d.cbf" % (i + 1)))
        writer = FullCBFWriter(imageset=imageset)
        frame_scan = Scan(
            (i + 1, i + 1),
            (10.0 + 0.25 * i, 0.25),
            exposure_times=flex.double([0.5]),
            epochs=flex.double([1.7e9 + i]),
        )
        cbf = writer.get_cbf_handle(
            index=0, header_only=True, goniometer=gonio, scan=frame_scan
        )
        writer.add_data_to_cbf(cbf, data=tuple(data))
        writer.write_cbf(filename, cbf=cbf)
        filenames.append(filename)

    experiments = ExperimentListFactory.from_filenames(filenames)
    assert len(experiments) == 1
    expt = experiments[0]
    assert expt.imageset.__class__.__name__ == "ImageSequence"
    assert expt.goniometer.is_similar_to(gonio)
    assert expt.scan.get_image_range() == scan.get_image_range()
    assert expt.scan.get_oscillation() == pytest.approx(scan.get_oscillation())
    assert list(expt.scan.get_exposure_times()) == list(scan.get_exposure_times())
    assert [p.get_image_size() for p in expt.detector] == sizes
    assert expt.detector.is_similar_to(det)
