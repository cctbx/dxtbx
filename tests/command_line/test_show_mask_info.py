from __future__ import annotations

from dxtbx.command_line import show_mask_info


def test_show_mask_info(dials_data, capsys):
    data = dials_data("image_examples") / "dectris_eiger_master.h5"

    show_mask_info.run([str(data)])
    captured = capsys.readouterr()
    assert not captured.err
    assert "Module 0 has 637992 masked pixels of 10166590" in captured.out
