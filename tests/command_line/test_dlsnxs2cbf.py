import procrunner
import pytest


@pytest.mark.xfail(
    "os.name == 'nt'", reason="crashes python process on Windows", run=False
)
def test_dlsnxs2cbf(dials_data, tmpdir):
    screen = dials_data("thaumatin_eiger_screen")
    master = screen.join("Therm_6_1_master.h5")
    result = procrunner.run(
        ["dxtbx.dlsnxs2cbf", master, "junk_%04d.cbf"], working_directory=tmpdir
    )
    assert not result.returncode and not result.stderr

    expected_output = "\n".join("junk_%04d.cbf" % j for j in (1, 2, 3))

    for record in expected_output.split("\n"):
        assert record.strip().encode("latin-1") in result.stdout, record

    # check files on disk

    for j in (1, 2, 3):
        assert tmpdir.join("junk_%04d.cbf" % j).check()
