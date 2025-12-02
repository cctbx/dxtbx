from __future__ import annotations

from dxtbx.command_line import print_header


def test_print_header(dials_data, capsys):
    screen = dials_data("thaumatin_eiger_screen")
    master = screen / "Therm_6_1_master.h5"
    print_header.run([str(master)])

    expected_output = [
        f"=== {master} ===",
        "Using header reader: FormatNXmxDLS16M",
    ]

    captured = capsys.readouterr()
    assert not captured.err
    for record in expected_output:
        assert record.strip() in captured.out, record
