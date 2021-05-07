import dxtbx.data.beamline_defs as dxbd


def test_lookup_unknown_detector():
    n = dxbd.get_beamline_definition(
        "This is a detector serial number that does not exist"
    )
    assert "Dummy CIF" in str(n), n
    assert str(n.CIF_block()) == ""
    assert n.CIF_block().__module__ == "iotbx.cif.model"
    assert str(n.mmCIF_block()) == ""
    assert n.mmCIF_block().__module__ == "iotbx.cif.model"


def test_lookup_known_detector():
    n = dxbd.get_beamline_definition("PILATUS 2M, S/N 24-0107 Diamond")
    assert "Dummy" not in str(n), n
    assert str(n) != ""
    cif = n.CIF_block()
    assert cif.__module__ == "iotbx.cif.model"
    assert "_diffrn_source" in str(cif)
    mmcif = n.mmCIF_block()
    assert mmcif.__module__ == "iotbx.cif.model"
    assert "_diffrn_source.source" in str(mmcif)
