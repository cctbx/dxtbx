import h5py

from dxtbx.format import nexus


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
