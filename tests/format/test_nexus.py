import h5py
import pytest

from dxtbx.format import nexus


def test_data_factory_nxs(dials_data):
    nxs_file = dials_data("vmxi_thaumatin") / "image_15799.nxs"
    with h5py.File(nxs_file, "r") as fh:
        data = fh["/entry/data"]
        data_factory = nexus.DataFactory(nexus.NXdata(data))
        for i in (0, data["data"].shape[0] - 1):
            assert data_factory[i].size()
        # If we try and access an image beyond the size of the dataset it should fail.
        # See https://github.com/cctbx/dxtbx/issues/285
        with pytest.raises(IndexError):
            data_factory[data["data"].shape[0]]


def test_find_entries(tmp_path):
    h5_file = tmp_path / "tmp.nxs"
    names = ("/entry", "/another_entry", "/yet_another_entry")
    with h5py.File(h5_file, "w") as fw:
        for name in names:
            group = fw.create_group(name)
            group.attrs["NX_class"] = "NXentry"
            group["definition"] = "NXmx"

    with h5py.File(h5_file, "r") as fr:
        entries = nexus.find_entries(fr)
        assert len(entries) == len(names)
        assert set(entry.name for entry in entries) == set(names)


def test_find_entries_no_entry(tmp_path):
    h5_file = tmp_path / "tmp.nxs"
    with h5py.File(h5_file, "w") as fw:
        # NX_class == NXdata
        group = fw.create_dataset("/data", (100,))
        group.attrs["NX_class"] = "NXdata"
        # NX_class == NXentry but no definition provided - skip this entry
        group = fw.create_group("/no_definition")
        group.attrs["NX_class"] = "NXentry"
        # NX_class == NXentry but definition != NXmx
        group = fw.create_group("/wrong_definition")
        group.attrs["NX_class"] = "NXentry"
        group["definition"] = "NXmagic"

    with h5py.File(h5_file, "r") as fr:
        entries = nexus.find_entries(fr)
        assert len(entries) == 0


def test_find_entries_empty_file(tmp_path):
    h5_file = tmp_path / "tmp.nxs"
    with h5py.File(h5_file, "w") as fw:
        fw.close()

    with h5py.File(h5_file, "r") as fr:
        entries = nexus.find_entries(fr)
        assert len(entries) == 0


def test_find_class_softlink(tmp_path):
    # 1) Create Nxbeam in /entry/sample/beam
    # 2) Store softlink to this beam in /entry/instrument/beam
    # 3) Call find_class on /entry/instrument and assert that it correctly follows the
    #    softlink
    h5_file = tmp_path / "tmp.nxs"
    with h5py.File(h5_file, "w") as fw:
        entry = fw.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        entry["definition"] = "NXmx"
        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        beam = sample.create_group("beam")
        beam.attrs["NX_class"] = "NXbeam"
        instrument = entry.create_group("instrument")
        instrument.attrs["NX_class"] = "NXinstrument"
        instrument["beam"] = h5py.SoftLink(beam.name)
        # Ensure that find_class gracefully handles broken soft/external links
        instrument["broken_link"] = h5py.SoftLink("/entry/nonsense")

    with h5py.File(h5_file, "r") as fr:
        assert len(nexus.find_class(fr["/entry/instrument"], "NXbeam")) == 1
        assert len(nexus.find_class(fr["/entry/sample"], "NXbeam")) == 1
