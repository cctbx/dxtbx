"""Tests for dxbtx.format.FormatNXmx format classes."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from dxtbx.format.FormatNXmx import FormatNXmx, FormatNXmxI19_2
from dxtbx.model.experiment_list import ExperimentListFactory


@pytest.fixture
def nxmx_example_on_disk(tmp_path, nxmx_example):
    """
    Copy the in-memory example NXmx data to a file on disk.

    Since ExperimentListFactory.from_filenames can only read from disk, and requires
    the presence of a pixel mask, create a modified on-disk copy of the NXmx example
    data.
    """
    with nxmx_example as f, h5py.File(tmp_path / "test.h5", "w") as g:
        f.copy("entry", g)

        # Define a pixel mask, as required for creating a DXTBX model.
        pixel_mask = "/entry/instrument/detector/pixel_mask"
        g.create_dataset(pixel_mask, data=np.zeros((4400, 4200)))

        yield g


parameters = [("DIAMOND BEAMLINE I19-2", FormatNXmxI19_2), ("DIAD", FormatNXmx)]


@pytest.mark.parametrize("instrument, format_class", parameters, ids=["I19-2", "DIAD"])
def test_Format_NXmx(nxmx_example_on_disk, instrument, format_class):
    """Check the right format class is used for the specified instrument."""
    name = "/entry/instrument/name"
    with nxmx_example_on_disk as g:
        del g[name]
        g.create_dataset(name, data=instrument)
        (experiment,) = ExperimentListFactory.from_filenames([g.filename])
        assert experiment.imageset.get_format_class() == format_class
