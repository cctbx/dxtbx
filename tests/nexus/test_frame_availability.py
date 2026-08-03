"""Tests for the live frame-availability helpers in dxtbx.nexus.

These model the way DLS NXmx masters are written during a live collection: the master
holds a virtual dataset (VDS) declaring the *full* planned frame count up front, whose
sources are external ``data_*.h5`` files. Crucially, a detector allocates each source
dataset at its full block extent from creation (``shape[0]`` is the planned length before
any pixel is written) and fills it one chunk per frame via SWMR. Availability must
therefore be judged from written chunks, not the (static) dataset extent.
"""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from dxtbx.nexus import get_available_frame_count, get_frame_counts

FRAME = (4, 4)


def _build_master(tmp_path, block_lens, stem="sim"):
    """Create an NXmx-like master with a '.'-self + ExternalLink VDS over preallocated,
    fixed-extent block files (one chunk per frame), mimicking a real detector layout.

    Returns (master_path, block_files) where block_files is an ordered list of
    (path, dset_name, block_len).
    """
    master = tmp_path / f"{stem}.nxs"
    total = sum(block_lens)
    blocks = []
    with h5py.File(master, "w", libver="latest") as m:
        entry = m.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        entry.create_dataset("definition", data="NXmx")
        data = entry.create_group("data")
        data.attrs["NX_class"] = "NXdata"
        data.attrs["signal"] = "data"

        layout = h5py.VirtualLayout(shape=(total,) + FRAME, dtype="uint32")
        start = 0
        for i, blen in enumerate(block_lens, 1):
            name = f"{stem}_{i:06d}.h5"
            path = tmp_path / name
            with h5py.File(path, "w", libver="latest") as sf:
                # Full extent up front, chunk-per-frame -- no chunks written yet.
                sf.create_dataset(
                    "data",
                    shape=(blen,) + FRAME,
                    chunks=(1,) + FRAME,
                    dtype="uint32",
                )
            link_name = f"data_{i:06d}"
            data[link_name] = h5py.ExternalLink(name, "data")
            layout[start : start + blen] = h5py.VirtualSource(
                ".", f"entry/data/{link_name}", shape=(blen,) + FRAME
            )
            blocks.append((path, "data", blen))
            start += blen
        data.create_virtual_dataset("data", layout, fillvalue=0)
    return master, blocks


def _write_frames(path, dset_name, n):
    """Write n leading frame chunks into a preallocated block file via SWMR."""
    with h5py.File(path, "r+", libver="latest") as f:
        f.swmr_mode = True
        dset = f[dset_name]
        for i in range(n):
            dset[i] = np.zeros(FRAME, dtype="uint32")
            dset.flush()


def test_preallocated_extent_reports_zero_until_written(tmp_path):
    """A freshly-created, full-extent block file must report 0 available, not its extent.

    This is the regression: reading ``shape[0]`` reported every block as fully available
    the instant its file existed, before any frame was written.
    """
    master, _ = _build_master(tmp_path, [5, 5])
    assert get_frame_counts(master) == (0, 10)
    assert get_available_frame_count(master) == 0


def test_partial_and_contiguous_counting(tmp_path):
    master, blocks = _build_master(tmp_path, [5, 5])
    # Partway through the first block.
    _write_frames(blocks[0][0], blocks[0][1], 3)
    assert get_frame_counts(master) == (3, 10)
    # First block complete, but nothing in the second yet.
    _write_frames(blocks[0][0], blocks[0][1], 5)
    assert get_frame_counts(master) == (5, 10)
    # Into the second block.
    _write_frames(blocks[1][0], blocks[1][1], 2)
    assert get_frame_counts(master) == (7, 10)
    # Fully collected.
    _write_frames(blocks[1][0], blocks[1][1], 5)
    assert get_frame_counts(master) == (10, 10)


def test_missing_source_file_stops_contiguous_count(tmp_path):
    """A block file that has not been created yet (rolling-window pre-creation) caps the
    contiguous available count without raising."""
    master, blocks = _build_master(tmp_path, [5, 5])
    _write_frames(blocks[0][0], blocks[0][1], 5)
    blocks[1][0].unlink()  # second block not yet on disk
    assert get_frame_counts(master) == (5, 10)


def test_file_close_method(tmp_path):
    """With method='file_close' a source counts only once it can be opened without SWMR."""
    master, blocks = _build_master(tmp_path, [5, 5])
    _write_frames(blocks[0][0], blocks[0][1], 5)
    # Files are closed after writing, so the first (full) block is available.
    assert get_available_frame_count(master, method="file_close") == 5


def test_non_virtual_dataset_is_fully_available(tmp_path):
    master = tmp_path / "plain.nxs"
    with h5py.File(master, "w") as m:
        entry = m.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        entry.create_dataset("definition", data="NXmx")
        data = entry.create_group("data")
        data.attrs["NX_class"] = "NXdata"
        data.attrs["signal"] = "data"
        data.create_dataset("data", data=np.zeros((7,) + FRAME, dtype="uint32"))
    assert get_frame_counts(master) == (7, 7)


def test_unreadable_master_returns_zeroes(tmp_path):
    assert get_frame_counts(tmp_path / "does_not_exist.nxs") == (0, 0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
