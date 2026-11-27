from __future__ import annotations

import h5py
import numpy as np
import pytest

from dxtbx.format.FormatESSNMX import get_tof_bin_edges

N_SLICES = 5
EDGES_NS = 73_588_426.0 + 300_000.0 * np.arange(N_SLICES + 1)
CENTRES_NS = 0.5 * (EDGES_NS[:-1] + EDGES_NS[1:])


def make_panel(tmp_path, **datasets):
    f = h5py.File(tmp_path / "panel.h5", "w")
    panel = f.create_group("entry/instrument/detector_panel_0")
    panel.create_dataset("data", data=np.zeros((4, 4, N_SLICES), dtype=np.int32))
    for name, values in datasets.items():
        panel.create_dataset(name, data=values)
    return f, panel


def test_edges_in_time_of_flight(tmp_path):
    # Older essnmx output: time_of_flight holds the n + 1 edges
    f, panel = make_panel(tmp_path, time_of_flight=EDGES_NS)
    np.testing.assert_allclose(get_tof_bin_edges(panel, 1e-3), EDGES_NS * 1e-3)
    f.close()


def test_centres_in_time_of_flight_with_edges(tmp_path):
    # Current essnmx output: n centres in time_of_flight, edges alongside
    f, panel = make_panel(
        tmp_path, time_of_flight=CENTRES_NS, original_time_edges=EDGES_NS
    )
    edges = get_tof_bin_edges(panel, 1e-3)
    assert edges.size == N_SLICES + 1
    np.testing.assert_allclose(edges, EDGES_NS * 1e-3)
    np.testing.assert_allclose(0.5 * (edges[:-1] + edges[1:]), CENTRES_NS * 1e-3)
    f.close()


def test_centres_only(tmp_path):
    f, panel = make_panel(tmp_path, time_of_flight=CENTRES_NS)
    np.testing.assert_allclose(get_tof_bin_edges(panel, 1e-3), EDGES_NS * 1e-3)
    f.close()


def test_mismatched_lengths_raise(tmp_path):
    f, panel = make_panel(tmp_path, time_of_flight=EDGES_NS[:3])
    with pytest.raises(ValueError):
        get_tof_bin_edges(panel, 1e-3)
    f.close()
