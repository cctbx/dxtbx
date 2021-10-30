import procrunner


def test_show_mask_info(dials_data):

    data = dials_data("image_examples") / "dectris_eiger_master.h5"

    result = procrunner.run(["dxtbx.show_mask_info", data.strpath])
    assert not result.returncode and not result.stderr

    assert b"Module 0 has 637992 masked pixels of 10166590" in result.stdout
