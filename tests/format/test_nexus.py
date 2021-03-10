from pathlib import Path

import h5py
import pytest

from dxtbx.format import nexus

# Some names for sample NXentries in a dummy NeXus file.
entry_names = "/entry", "/another_entry", "/yet_another_entry"


@pytest.fixture(scope="session")
def nexus_file(tmp_path_factory) -> Path:
    """
    Create a dummy NXmx-flavoured NeXus HDF5 file to exercise parts of format/nexus.py.

    The file contains
        - Three NXentries.
        - Two NXsamples within the NXentry 'entry'.
        - A NXinstrument within 'entry'.
        - A NXbeam within the NXsample 'sample', with a soft link from the NXinstrument.
        - A broken soft link from the NXinstrument.
        - A HDF5 group and a data set without NXclass attributes, both within 'entry'.

    Args:
        tmp_path_factory:  Pytest session-scoped temporary directory fixture.

    Returns:
        Path to the dummy NeXus file.
    """
    h5_file = tmp_path_factory.mktemp("dxtbx_dummy_nexus", numbered=False) / "tmp.nxs"

    with h5py.File(h5_file, "w") as fw:
        # Create some top-level NXmx NXentries.
        for name in entry_names:
            group = fw.create_group(name)
            group.attrs["NX_class"] = "NXentry"
            group["definition"] = "NXmx"
        entry = fw["entry"]

        # Create two NXsamples in the NXentry entry.
        for name in "sample", "another_sample":
            dataset = entry.create_group(name)
            dataset.attrs["NX_class"] = "NXsample"
        sample = entry["sample"]

        # Create a NXinstrument in the NXentry entry.
        instrument = entry.create_group("instrument")
        instrument.attrs["NX_class"] = "NXinstrument"

        # Create a NXbeam in the NXsample sample.
        beam = sample.create_group("beam")
        beam.attrs["NX_class"] = "NXbeam"
        # Link to it from the NXinstrument instrument.
        instrument["beam"] = h5py.SoftLink(beam.name)

        # Ensure that find_classes gracefully handles broken soft/external links
        instrument["broken_link"] = h5py.SoftLink("/entry/nonsense")

        # Ensure that find_classes copes with non-NX_class groups and data sets.
        entry.create_group("group")
        entry.create_dataset("dataset", data="Test data set.")

    return h5_file


def test_data_factory_nxs(dials_data):
    nxs_file = dials_data("vmxi_thaumatin") / "image_15799.nxs"
    with h5py.File(nxs_file) as fh:
        data = fh["/entry/data"]
        data_factory = nexus.DataFactory(nexus.NXdata(data))
        for i in (0, data["data"].shape[0] - 1):
            assert data_factory[i].size()
        # If we try and access an image beyond the size of the dataset it should fail.
        # See https://github.com/cctbx/dxtbx/issues/285
        with pytest.raises(IndexError):
            data_factory[data["data"].shape[0]]


def test_find_entries(nexus_file):
    with h5py.File(nexus_file, "r") as fr:
        entries = nexus.find_entries(fr)
        assert len(entries) == len(entry_names)
        assert set(entry.name for entry in entries) == set(entry_names)


def test_find_entries_no_entry(tmp_path):
    h5_file = tmp_path / "tmp.nxs"
    with h5py.File(h5_file, "w") as fw:
        # NX_class == NXdata
        group = fw.create_dataset("data", (100,))
        group.attrs["NX_class"] = "NXdata"
        # NX_class == NXentry but no definition provided - skip this entry
        group = fw.create_group("no_definition")
        group.attrs["NX_class"] = "NXentry"
        # NX_class == NXentry but definition != NXmx
        group = fw.create_group("wrong_definition")
        group.attrs["NX_class"] = "NXentry"
        group["definition"] = "NXmagic"

    with h5py.File(h5_file, "r") as fr:
        entries = nexus.find_entries(fr)
        assert len(entries) == 0


def test_find_entries_empty_file(tmp_path):
    h5_file = tmp_path / "tmp.nxs"
    with h5py.File(h5_file, "w"):
        pass

    with h5py.File(h5_file, "r") as fr:
        entries = nexus.find_entries(fr)
        assert len(entries) == 0


def test_find_classes(nexus_file):
    """
    Test dxtbx.format.nexus.find_classes.

    1. Call find_classes on the main NXentry of the dummy NeXus file.
    2. Ensure it finds one instrument, two samples and three children that do not
       belong to a NX_class.
    3. Ensure that the children that don't belong to a NX_class are 'definition',
       'group' and 'dataset'.
    """
    with h5py.File(nexus_file, "r") as fr:
        samples, instruments, classless = nexus.find_classes(
            fr["entry"], "NXsample", "NXinstrument", None
        )
        assert len(instruments) == 1
        assert len(samples) == 2
        assert len(classless) == 3
        expected_classless = {"/entry/definition", "/entry/group", "/entry/dataset"}
        assert {value.name for value in classless} == expected_classless


def test_find_class(nexus_file):
    """
    Test dxtbx.format.nexus.find_class

    1. Ensure that nexus.find_class finds the single NXbeam child of /entry/sample.
    2. Ensure that the broken link /entry/instrument/broken_link exists but is
       ignored by nexus.find_class, even in a search for children without an NX_class.
    """
    with h5py.File(nexus_file, "r") as fr:
        beams = nexus.find_class(fr["entry/sample"], "NXbeam")
        assert len(beams) == 1
        assert beams[0].name == "/entry/sample/beam"

        instrument = fr["entry/instrument"]
        classless = nexus.find_class(instrument, None)
        assert "broken_link" in instrument and not classless


def test_find_class_softlink(nexus_file):
    """
    Test that nexus.find_class follows soft links.

    Call find_class on /entry/instrument and assert that it correctly follows the
    soft link to /entry/sample/beam.
    """
    with h5py.File(nexus_file, "r") as fr:
        assert len(nexus.find_class(fr["entry/instrument"], "NXbeam")) == 1
        assert len(nexus.find_class(fr["entry/sample"], "NXbeam")) == 1
