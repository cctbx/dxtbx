from __future__ import annotations

import collections
import errno
import os
import pickle
from unittest import mock

import pytest

from cctbx import sgtbx
from scitbx.array_family import flex

import dxtbx
import dxtbx.model.experiment_list
from dxtbx.format.Format import Format
from dxtbx.imageset import ImageSetFactory
from dxtbx.model import (
    Beam,
    BeamFactory,
    Crystal,
    Detector,
    Experiment,
    ExperimentList,
    ExperimentType,
    Goniometer,
    Scan,
    ScanFactory,
)

try:
    from ..dxtbx_model_ext import Probe
except ModuleNotFoundError:
    from dxtbx_model_ext import Probe  # type: ignore
from dxtbx.model.experiment_list import ExperimentListDict, ExperimentListFactory


@pytest.fixture(scope="session")
def centroid_test_data(dials_regression):
    return os.path.join(dials_regression, "centroid_test_data")


@pytest.fixture
def single_sequence_filenames(centroid_test_data):
    filenames = [
        os.path.join(centroid_test_data, f"centroid_000{i}.cbf") for i in range(1, 10)
    ]
    return filenames


@pytest.fixture
def multiple_sequence_filenames(centroid_test_data):
    filenames = [
        os.path.join(centroid_test_data, f"centroid_000{i}.cbf")
        for i in [1, 2, 3, 7, 8, 9]
    ]
    return filenames


@pytest.fixture
def all_image_examples(dials_regression):
    filenames = (
        ("ALS_1231", "q315r_lyso_1_001.img"),
        ("ALS_501", "als501_q4_1_001.img"),
        ("ALS_821", "q210_lyso_1_101.img"),
        ("ALS_831", "q315r_lyso_001.img"),
        ("APS_14BMC", "q315_1_001.img"),
        ("APS_17ID", "q210_1_001.img"),
        ("APS_19ID", "q315_unbinned_a.0001.img"),
        ("APS_22ID", "mar300.0001"),
        ("APS_23IDD", "mar300_1_E1.0001"),
        ("APS_24IDC", "pilatus_1_0001.cbf"),
        ("APS_24IDC", "q315_1_001.img"),
        ("CLS1_08ID1", "mar225_2_E0_0001.img"),
        ("DESY_ID141", "q210_2_001.img"),
        ("ESRF_BM14", "mar165_001.mccd"),
        ("ESRF_BM14", "mar225_1_001.mccd"),
        ("ESRF_ID231", "q315r_7_001.img"),
        ("RAXIS-HTC", "test1_lysozyme_0111060001.osc"),
        ("SLS_X06SA", "mar225_2_001.img"),
        ("SLS_X06SA", "pilatus6m_1_00001.cbf"),
        ("SRS_101", "mar225_001.img"),
        ("SRS_142", "q4_1_001.img"),
        ("SSRL_bl111", "mar325_1_001.mccd"),
        ("xia2", "merge2cbf_averaged_0001.cbf"),
    )
    return [os.path.join(dials_regression, "image_examples", *f) for f in filenames]


@pytest.fixture
def multiple_block_filenames(single_sequence_filenames, all_image_examples):
    return single_sequence_filenames + all_image_examples


def test_experiment_list_extend():
    """Check that the extend method of ExperimentList works."""
    # Create a minimal ExperimentList instance.
    expts = ExperimentList([Experiment()])
    # Try to extend it.
    expts.extend(expts)


def test_experiment_contains():
    # Create a load of models
    b1 = Beam()
    d1 = Detector()
    g1 = Goniometer()
    s1 = Scan()
    c1 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    # Create an experiment
    e = Experiment(
        beam=b1, detector=d1, goniometer=g1, scan=s1, crystal=c1, imageset=None
    )

    # Check experiment contains model
    assert b1 in e
    assert d1 in e
    assert g1 in e
    assert s1 in e
    assert c1 in e

    # Create a load of models that look the same but aren't
    b2 = Beam()
    d2 = Detector()
    g2 = Goniometer()
    s2 = Scan()
    c2 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    # Check experiment doesn't contain model
    assert b2 not in e
    assert d2 not in e
    assert g2 not in e
    assert s2 not in e
    assert c2 not in e


def test_experiment_equality():
    # Create a load of models
    b1 = Beam()
    d1 = Detector()
    g1 = Goniometer()
    s1 = Scan()
    c1 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    # Create a load of models that look the same but aren't
    b2 = Beam()
    d2 = Detector()
    g2 = Goniometer()
    s2 = Scan()
    c2 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    # Create a few experiments
    e1 = Experiment(
        beam=b1, detector=d1, goniometer=g1, scan=s1, crystal=c1, imageset=None
    )
    e2 = Experiment(
        beam=b1, detector=d1, goniometer=g1, scan=s1, crystal=c1, imageset=None
    )
    e3 = Experiment(
        beam=b2, detector=d2, goniometer=g2, scan=s2, crystal=c2, imageset=None
    )

    # Check e1 equals e2 and neither equals e3
    assert e1 == e2
    assert e1 != e3
    assert e2 != e3


def test_experiment_consistent(dials_data):
    # Create a sequence
    sequence_filenames = dials_data("centroid_test_data", pathlib=False).listdir(
        "centroid*.cbf"
    )
    sequence = ImageSetFactory.new(sorted(f.strpath for f in sequence_filenames))[0]

    # Create experiment with sequence and good scan
    e = Experiment(imageset=sequence, scan=sequence.get_scan())
    assert e.is_consistent()

    # Create experiment with sequence and defective scan
    scan = sequence.get_scan()
    scan.set_image_range((1, 1))
    e = Experiment(imageset=sequence, scan=scan)
    # assert not e.is_consistent()) # FIXME

    ## Create experiment with imageset and good scan
    # assert e.is_consistent()

    ## Create experiment with imageset and non-still scan
    # assert not e.is_consistent()

    ## Create experiment with imageset and scan with more than 1 image
    # assert not e.is_consistent()

    ## Create experiment with imageset and defective scan
    # assert not e.is_consistent()


def test_experimentlist_contains(experiment_list):
    # Check all the models are found
    for e in experiment_list:
        assert e.beam in experiment_list
        assert e.detector in experiment_list
        assert e.goniometer in experiment_list
        assert e.scan in experiment_list

    # Create some more models
    b = Beam()
    d = Detector()
    g = Goniometer()
    s = Scan()

    # Check that models not in are not found
    assert b not in experiment_list
    assert d not in experiment_list
    assert g not in experiment_list
    assert s not in experiment_list


def test_experimentlist_replace(experiment_list):
    # Get the models
    b = [e.beam for e in experiment_list]

    # Replace some models
    experiment_list.replace(b[0], b[1])
    assert experiment_list[0].beam is b[1]
    assert experiment_list[4].beam is b[1]

    # Replace again
    experiment_list[0].beam = b[0]
    experiment_list[4].beam = b[4]


def test_experimentlist_indices(experiment_list):
    # Get the models
    b = [e.beam for e in experiment_list]
    d = [e.detector for e in experiment_list]
    g = [e.goniometer for e in experiment_list]
    s = [e.scan for e in experiment_list]

    # Check indices of beams
    assert list(experiment_list.indices(b[0])) == [0, 4]
    assert list(experiment_list.indices(b[1])) == [1, 3]
    assert list(experiment_list.indices(b[2])) == [2]
    assert list(experiment_list.indices(b[3])) == [1, 3]
    assert list(experiment_list.indices(b[4])) == [0, 4]

    # Check indices of detectors
    assert list(experiment_list.indices(d[0])) == [0, 4]
    assert list(experiment_list.indices(d[1])) == [1, 3]
    assert list(experiment_list.indices(d[2])) == [2]
    assert list(experiment_list.indices(d[3])) == [1, 3]
    assert list(experiment_list.indices(d[4])) == [0, 4]

    # Check indices of goniometer
    assert list(experiment_list.indices(g[0])) == [0, 4]
    assert list(experiment_list.indices(g[1])) == [1, 3]
    assert list(experiment_list.indices(g[2])) == [2]
    assert list(experiment_list.indices(g[3])) == [1, 3]
    assert list(experiment_list.indices(g[4])) == [0, 4]

    # Check indices of scans
    assert list(experiment_list.indices(s[0])) == [0, 4]
    assert list(experiment_list.indices(s[1])) == [1, 3]
    assert list(experiment_list.indices(s[2])) == [2]
    assert list(experiment_list.indices(s[3])) == [1, 3]
    assert list(experiment_list.indices(s[4])) == [0, 4]

    # Check some models not in the list
    assert len(experiment_list.indices(Beam())) == 0
    assert len(experiment_list.indices(Detector())) == 0
    assert len(experiment_list.indices(Goniometer())) == 0
    assert len(experiment_list.indices(Scan())) == 0


def test_experimentlist_models(experiment_list):
    # Get all the unique models
    b = experiment_list.beams()
    d = experiment_list.detectors()
    g = experiment_list.goniometers()
    s = experiment_list.scans()

    # Check we have the expected number
    assert len(b) == 3
    assert len(d) == 3
    assert len(g) == 3
    assert len(s) == 3

    # Check we have the expected order
    assert b[0] == experiment_list[0].beam
    assert b[1] == experiment_list[1].beam
    assert b[2] == experiment_list[2].beam

    assert d[0] == experiment_list[0].detector
    assert d[1] == experiment_list[1].detector
    assert d[2] == experiment_list[2].detector

    assert g[0] == experiment_list[0].goniometer
    assert g[0] == experiment_list[0].goniometer
    assert g[1] == experiment_list[1].goniometer

    assert s[2] == experiment_list[2].scan
    assert s[1] == experiment_list[1].scan
    assert s[2] == experiment_list[2].scan


def test_experimentlist_to_dict(experiment_list):
    # Convert the list to a dictionary
    obj = experiment_list.to_dict()

    # Check this is the right object
    assert obj["__id__"] == "ExperimentList"

    # Check length of items
    assert len(obj["experiment"]) == 5
    assert len(obj["beam"]) == 3
    assert len(obj["detector"]) == 3
    assert len(obj["goniometer"]) == 3
    assert len(obj["scan"]) == 3

    # The expected models
    b = [0, 1, 2, 1, 0]
    d = [0, 1, 2, 1, 0]
    g = [0, 1, 2, 1, 0]
    s = [0, 1, 2, 1, 0]

    # Check all the experiments
    for i, eobj in enumerate(obj["experiment"]):
        assert eobj["__id__"] == "Experiment"
        assert eobj["beam"] == b[i]
        assert eobj["detector"] == d[i]
        assert eobj["goniometer"] == g[i]
        assert eobj["scan"] == s[i]


def test_experimentlist_where(experiment_list):
    for beam in experiment_list.beams():
        assert beam is not None
        for i in experiment_list.where(beam=beam):
            assert experiment_list[i].beam is beam
    for goniometer in experiment_list.goniometers():
        assert goniometer is not None
        for i in experiment_list.where(goniometer=goniometer):
            assert experiment_list[i].goniometer is goniometer
    for scan in experiment_list.scans():
        assert scan is not None
        for i in experiment_list.where(scan=scan):
            assert experiment_list[i].scan is scan
    for detector in experiment_list.detectors():
        assert detector is not None
        for i in experiment_list.where(detector=detector):
            assert experiment_list[i].detector is detector


@pytest.fixture
def experiment_list():
    # Initialise a list of experiments
    experiments = ExperimentList()

    # Create a few beams
    b1 = Beam()
    b2 = Beam()
    b3 = Beam()

    # Create a few detectors
    d1 = Detector()
    d2 = Detector()
    d3 = Detector()

    # Create a few goniometers
    g1 = Goniometer()
    g2 = Goniometer()
    g3 = Goniometer()

    # Create a few scans
    s1 = Scan()
    s2 = Scan()
    s3 = Scan()

    # Create a list of models
    b = [b1, b2, b3, b2, b1]
    d = [d1, d2, d3, d2, d1]
    g = [g1, g2, g3, g2, g1]
    s = [s1, s2, s3, s2, s1]
    ident = ["sausage", "eggs", "bacon", "toast", "beans"]

    # Populate with various experiments
    for i in range(5):
        experiments.append(
            Experiment(
                beam=b[i],
                detector=d[i],
                goniometer=g[i],
                scan=s[i],
                identifier=ident[i],
            )
        )

    return experiments


def test_experimentlist_factory_from_json(monkeypatch, dials_regression):

    # Get all the filenames
    filename1 = os.path.join(
        dials_regression, "experiment_test_data", "experiment_1.json"
    )
    filename3 = os.path.join(
        dials_regression, "experiment_test_data", "experiment_3.json"
    )
    filename4 = os.path.join(
        dials_regression, "experiment_test_data", "experiment_4.json"
    )

    # Read all the experiment lists in
    with monkeypatch.context() as m:
        m.setenv("DIALS_REGRESSION", dials_regression)
        el1 = ExperimentListFactory.from_json_file(filename1)
        el3 = ExperimentListFactory.from_json_file(filename3)
        el4 = ExperimentListFactory.from_json_file(filename4)

    # All the experiment lists should be the same length
    assert len(el1) == 1
    assert len(el1) == len(el3)
    assert len(el1) == len(el4)

    # Check all the models are the same
    for e in zip(el1, el3, el4):
        e1 = e[0]
        assert e1.imageset
        assert e1.beam
        assert e1.detector
        assert e1.goniometer
        assert e1.scan
        assert e1.crystal
        for ee in e[1:]:
            assert e1.imageset == ee.imageset
            assert e1.beam == ee.beam
            assert e1.detector == ee.detector
            assert e1.goniometer == ee.goniometer
            assert e1.scan == ee.scan
            assert e1.crystal == ee.crystal


def test_experimentlist_factory_from_pickle(monkeypatch, dials_regression):
    # Get all the filenames
    filename1 = os.path.join(
        dials_regression, "experiment_test_data", "experiment_1.json"
    )

    # Read all the experiment lists in
    with monkeypatch.context() as m:
        m.setenv("DIALS_REGRESSION", dials_regression)
        el1 = ExperimentListFactory.from_json_file(filename1)

    # Pickle then load again
    el2 = pickle.loads(pickle.dumps(el1))

    # All the experiment lists should be the same length
    assert len(el1) == 1
    assert len(el1) == len(el2)

    # Check all the models are the same
    for e1, e2 in zip(el1, el2):
        assert e1.imageset and e1.imageset == e2.imageset
        assert e1.beam and e1.beam == e2.beam
        assert e1.detector and e1.detector == e2.detector
        assert e1.goniometer and e1.goniometer == e2.goniometer
        assert e1.scan and e1.scan == e2.scan
        assert e1.crystal and e1.crystal == e2.crystal


def test_experimentlist_factory_from_args(monkeypatch, dials_regression):
    pytest.importorskip("dials")

    # Get all the filenames
    filenames = [
        os.path.join(dials_regression, "experiment_test_data", "experiment_1.json"),
        # os.path.join(dials_regression, 'experiment_test_data', 'experiment_2.json'),
        os.path.join(dials_regression, "experiment_test_data", "experiment_3.json"),
        os.path.join(dials_regression, "experiment_test_data", "experiment_4.json"),
    ]

    # Get the experiments from a list of filenames
    with monkeypatch.context() as m:
        m.setenv("DIALS_REGRESSION", dials_regression)
        experiments = ExperimentListFactory.from_args(filenames)

    assert len(experiments) == 3
    for experiment in experiments:
        assert experiment.imageset
        assert experiment.beam
        assert experiment.detector
        assert experiment.goniometer
        assert experiment.scan


def test_experimentlist_factory_from_imageset():
    imageset = Format.get_imageset(["filename.cbf"], as_imageset=True)
    imageset.set_beam(Beam(), 0)
    imageset.set_detector(Detector(), 0)

    crystal = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    experiments = ExperimentListFactory.from_imageset_and_crystal(imageset, crystal)

    assert len(experiments) == 1
    assert experiments[0].imageset
    assert experiments[0].beam
    assert experiments[0].detector is not None
    assert experiments[0].crystal


def test_experimentlist_factory_from_sequence():
    filenames = ["filename_%01d.cbf" % (i + 1) for i in range(0, 2)]

    imageset = Format.get_imageset(
        filenames,
        beam=Beam(),
        detector=Detector(),
        goniometer=Goniometer(),
        scan=Scan((1, 2), (0, 1)),
        as_sequence=True,
    )

    crystal = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    experiments = ExperimentListFactory.from_imageset_and_crystal(imageset, crystal)

    assert len(experiments) == 1
    assert experiments[0].imageset
    assert experiments[0].beam
    assert experiments[0].detector is not None
    assert experiments[0].goniometer
    assert experiments[0].scan
    assert experiments[0].crystal


def test_experimentlist_dumper_dump_formats(monkeypatch, dials_regression, tmp_path):
    # Get all the filenames
    filename1 = os.path.join(
        dials_regression, "experiment_test_data", "experiment_1.json"
    )

    # Read all the experiment lists in
    with monkeypatch.context() as m:
        m.setenv("DIALS_REGRESSION", dials_regression)
        elist1 = ExperimentListFactory.from_json_file(filename1)

    # Dump as JSON file and reload
    filename = tmp_path / "temp1.json"
    elist1.as_json(filename)
    elist2 = ExperimentListFactory.from_json_file(filename)
    check(elist1, elist2)

    # Dump as split JSON file and reload
    filename = tmp_path / "temp2.json"
    elist1.as_json(filename, split=True)
    elist2 = ExperimentListFactory.from_json_file(filename)
    check(elist1, elist2)


def test_experimentlist_dumper_dump_scan_varying(
    monkeypatch, dials_regression, tmp_path
):
    # Get all the filenames
    filename1 = os.path.join(
        dials_regression, "experiment_test_data", "experiment_1.json"
    )

    # Read the experiment list in
    with monkeypatch.context() as m:
        m.setenv("DIALS_REGRESSION", dials_regression)
        elist1 = ExperimentListFactory.from_json_file(filename1)

    # Make trivial scan-varying models
    crystal = elist1[0].crystal
    beam = elist1[0].beam
    goniometer = elist1[0].goniometer
    crystal.set_A_at_scan_points([crystal.get_A()] * 5)

    cov_B = flex.double([1e-5] * 9 * 9)
    crystal.set_B_covariance(cov_B)
    cov_B.reshape(flex.grid(1, 9, 9))
    cov_B_array = flex.double(flex.grid(5, 9, 9))
    for i in range(5):
        cov_B_array[i : (i + 1), :, :] = cov_B
    crystal.set_B_covariance_at_scan_points(cov_B_array)

    beam.set_s0_at_scan_points([beam.get_s0()] * 5)
    goniometer.set_setting_rotation_at_scan_points(
        [goniometer.get_setting_rotation()] * 5
    )

    # Dump as JSON file and reload
    filename = tmp_path / "temp.json"
    elist1.as_json(filename)
    elist2 = ExperimentListFactory.from_json_file(filename)
    check(elist1, elist2)


def test_experimentlist_dumper_dump_empty_sequence(tmp_path):
    filenames = [tmp_path / f"filename_{i}.cbf" for i in range(1, 3)]

    imageset = Format.get_imageset(
        filenames,
        beam=Beam((1, 0, 0)),
        detector=Detector(),
        goniometer=Goniometer(),
        scan=Scan((1, 2), (0.0, 1.0)),
        as_sequence=True,
    )

    crystal = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    experiments = ExperimentListFactory.from_imageset_and_crystal(imageset, crystal)

    filename = tmp_path / "temp.json"
    experiments.as_json(filename)
    experiments2 = ExperimentListFactory.from_json_file(filename, check_format=False)
    check(experiments, experiments2)


def test_experimentlist_dumper_dump_with_lookup(dials_regression, tmp_path):
    filename = os.path.join(
        dials_regression, "centroid_test_data", "experiments_with_lookup.json"
    )

    experiments = ExperimentListFactory.from_json_file(filename, check_format=True)

    imageset = experiments[0].imageset
    assert not imageset.external_lookup.mask.data.empty()
    assert not imageset.external_lookup.gain.data.empty()
    assert not imageset.external_lookup.pedestal.data.empty()
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.tile(0).data().all_eq(True)
    assert imageset.external_lookup.gain.data.tile(0).data().all_eq(1)
    assert imageset.external_lookup.pedestal.data.tile(0).data().all_eq(0)

    filename = tmp_path / "temp.json"
    experiments.as_json(filename)

    experiments = ExperimentListFactory.from_json_file(filename, check_format=True)

    imageset = experiments[0].imageset
    assert not imageset.external_lookup.mask.data.empty()
    assert not imageset.external_lookup.gain.data.empty()
    assert not imageset.external_lookup.pedestal.data.empty()
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.tile(0).data().all_eq(True)
    assert imageset.external_lookup.gain.data.tile(0).data().all_eq(1)
    assert imageset.external_lookup.pedestal.data.tile(0).data().all_eq(0)


def test_experimentlist_dumper_dump_with_bad_lookup(dials_data, tmpdir):
    filename = (
        dials_data("centroid_test_data", pathlib=True)
        / "experiments_with_bad_lookup.json"
    )
    experiments = ExperimentListFactory.from_json_file(filename, check_format=False)

    imageset = experiments[0].imageset
    assert imageset.external_lookup.mask.data.empty()
    assert imageset.external_lookup.gain.data.empty()
    assert imageset.external_lookup.pedestal.data.empty()
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None

    filename = tmpdir / "temp.json"
    experiments.as_json(filename)

    experiments = ExperimentListFactory.from_json_file(filename, check_format=False)

    imageset = experiments[0].imageset
    assert imageset.external_lookup.mask.data.empty()
    assert imageset.external_lookup.gain.data.empty()
    assert imageset.external_lookup.pedestal.data.empty()
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None


def test_experimentlist_with_identifiers():
    # Initialise a list of experiments
    experiments = ExperimentList()

    experiments.append(
        Experiment(beam=Beam(s0=(0, 0, -1)), detector=Detector(), identifier="bacon")
    )

    experiments.append(
        Experiment(beam=Beam(s0=(0, 0, -1)), detector=Detector(), identifier="sausage")
    )

    with pytest.raises(Exception):
        experiments.append(
            Experiment(beam=Beam(), detector=Detector(), identifier="bacon")
        )

    d = experiments.to_dict()
    e2 = ExperimentListDict(d).decode()

    assert experiments[0].identifier == e2[0].identifier
    assert experiments[1].identifier == e2[1].identifier

    assert tuple(experiments.identifiers()) == ("bacon", "sausage")
    experiments[0].identifier = "spam"
    assert tuple(experiments.identifiers()) == ("spam", "sausage")

    experiments.append(Experiment(identifier="bacon"))
    experiments.select_on_experiment_identifiers(["spam", "bacon"])
    assert list(experiments.identifiers()) == ["spam", "bacon"]
    experiments.append(Experiment(identifier="ham"))
    experiments.append(Experiment(identifier="jam"))
    experiments.remove_on_experiment_identifiers(["spam", "jam"])
    assert list(experiments.identifiers()) == ["bacon", "ham"]


def test_load_models(dials_data):
    pytest.importorskip("h5py")
    filename = (
        dials_data("image_examples", pathlib=True) / "SACLA-MPCCD-run266702-0-subset.h5"
    )

    # Test different ways of loading the data
    waves1, waves2, waves3, waves4, waves5 = [], [], [], [], []
    oris1, oris2, oris3, oris4, oris5 = [], [], [], [], []

    img = dxtbx.load(filename)

    # Test using dxtbx directly
    for i in range(img.get_num_images()):
        waves1.append(img.get_beam(i).get_wavelength())
        oris1.append(img.get_detector(i)[0].get_origin())

    # Test using the imageset clases
    imageset = img.get_imageset(str(filename))
    for i in range(len(imageset)):
        waves2.append(imageset.get_beam(i).get_wavelength())
        oris2.append(imageset.get_detector(i)[0].get_origin())

    # Test using imageset subsets
    imageset = img.get_imageset(str(filename))
    for i in range(len(imageset)):
        subset = imageset[i : i + 1]
        waves3.append(subset.get_beam(0).get_wavelength())
        oris3.append(subset.get_detector(0)[0].get_origin())

    # Test using pre-loaded experiments
    experiments = ExperimentListFactory.from_filenames([filename])
    for experiment in experiments:
        waves4.append(experiment.beam.get_wavelength())
        oris4.append(experiment.detector[0].get_origin())

    # Test using post-loaded experiments
    experiments = ExperimentListFactory.from_filenames([filename], load_models=False)
    for experiment in experiments:
        assert experiment.beam is None and experiment.detector is None
        experiment.load_models()
        waves5.append(experiment.beam.get_wavelength())
        oris5.append(experiment.detector[0].get_origin())

    for w1, w2, w3, w4, w5 in zip(waves1, waves2, waves3, waves4, waves5):
        assert w1 == w2 == w3 == w4 == w5

    for o1, o2, o3, o4, o5 in zip(oris1, oris2, oris3, oris4, oris5):
        assert o1 == o2 == o3 == o4 == o5


def test_partial_missing_model_serialization():
    goniometer = Goniometer()
    elist = ExperimentList([Experiment(), Experiment(goniometer=goniometer)])
    elist_ = ExperimentListFactory.from_dict(elist.to_dict())
    check(elist, elist_)


def test_experiment_type():
    experiment = Experiment()
    assert experiment.get_type() == ExperimentType.STILL
    experiment.goniometer = Goniometer()
    assert experiment.get_type() == ExperimentType.STILL
    experiment.scan = Scan()
    assert experiment.get_type() == ExperimentType.STILL
    experiment.scan = Scan((1, 1000), (0, 0.05))
    assert experiment.get_type() == ExperimentType.ROTATION
    # Specifically test the bug from dxtbx#4 triggered by ending on 0°
    experiment.scan = Scan((1, 1800), (-90, 0.05))
    assert experiment.get_type() == ExperimentType.ROTATION

    experiment.beam = BeamFactory.make_polychromatic_beam(
        direction=(0, 0, -1),
        sample_to_source_distance=(100),
        probe=Probe.xray,
        wavelength_range=(1, 10),
    )

    assert experiment.get_type() == ExperimentType.LAUE

    experiment.scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"time_of_flight": list(range(10))}
    )
    assert experiment.get_type() == ExperimentType.TOF

    experiment.beam = Beam()
    experiment.scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"other_property": list(range(10))}
    )
    assert experiment.get_type() == ExperimentType.STILL


def check(el1, el2):
    # All the experiment lists should be the same length
    assert len(el1) == len(el2)

    # Check all the models are the same
    for e1, e2 in zip(el1, el2):
        assert compare_experiment(e1, e2)


def compare_experiment(exp1, exp2):
    return (
        exp1.imageset == exp2.imageset
        and exp1.beam == exp2.beam
        and exp1.detector == exp2.detector
        and exp1.goniometer == exp2.goniometer
        and exp1.scan == exp2.scan
        and exp1.profile == exp2.profile
        and exp1.scaling_model == exp2.scaling_model
        and exp1.identifier == exp2.identifier
    )


def test_experimentlist_from_file(monkeypatch, dials_regression, tmpdir):
    # With the default check_format=True this file should fail to load with an
    # appropriate error as we can't find the images on disk
    with monkeypatch.context() as m:
        m.delenv("DIALS_REGRESSION", raising=False)
        with pytest.raises(IOError) as e:
            exp_list = ExperimentList.from_file(
                os.path.join(
                    dials_regression, "experiment_test_data", "experiment_1.json"
                )
            )
    assert e.value.errno == errno.ENOENT
    assert "No such file or directory" in str(e.value)
    assert "centroid_0001.cbf" in str(e.value)

    # Setting check_format=False should allow the file to load
    exp_list = ExperimentList.from_file(
        os.path.join(dials_regression, "experiment_test_data", "experiment_1.json"),
        check_format=False,
    )
    assert len(exp_list) == 1
    assert exp_list[0].beam

    # This allows expansion of environment variables in regression files, enabling the
    # file to load with check_format=True

    with monkeypatch.context() as m:
        m.setenv("DIALS_REGRESSION", dials_regression)
        exp_list = ExperimentList.from_file(
            os.path.join(dials_regression, "experiment_test_data", "experiment_1.json")
        )
    assert len(exp_list) == 1
    assert exp_list[0].beam


def test_experimentlist_imagesequence_stills(dials_data):
    filenames = [
        str(
            dials_data("thaumatin_grid_scan", pathlib=True)
            / f"thau_3_2_{i:04d}.cbf.bz2"
        )
        for i in range(1, 4)
    ]
    experiments = ExperimentListFactory.from_filenames(filenames)

    assert len(experiments) == 3
    assert len(experiments.imagesets()) == 1

    # Convert experiment list to dict
    d = experiments.to_dict()

    # Decode the dict to get a new experiment list
    experiments2 = ExperimentListDict(d).decode()

    # Verify that this experiment is as we expect
    assert len(experiments2) == 3
    assert len(experiments2.imagesets()) == 1
    assert len(experiments2.goniometers()) == 1
    assert len(experiments2.detectors()) == 1
    assert len(experiments2.beams()) == 1
    assert len(experiments2.scans()) == 3
    for expt in experiments2:
        assert expt.imageset is experiments2.imagesets()[0]


def test_experimentlist_imagesequence_decode(mocker):
    # These models are shared between experiments
    beam = Beam(s0=(0, 0, -1))
    detector = Detector()
    gonio = Goniometer()

    # Construct the experiment list
    experiments = ExperimentList()
    for i in range(3):
        experiments.append(
            Experiment(
                beam=beam,
                detector=detector,
                scan=ScanFactory.make_scan(
                    image_range=(i + 1, i + 1),
                    exposure_times=[1],
                    oscillation=(0, 0),
                    epochs=[0],
                ),
                goniometer=gonio,
            )
        )

    # Convert experiment list to dict and manually insert a shared imageset
    d = experiments.to_dict()
    d["imageset"].append(
        {"__id__": "ImageSequence", "template": "Puck3_10_1_####.cbf.gz"}
    )
    for e in d["experiment"]:
        e["imageset"] = 0

    # Monkeypatch this function as we don't actually have an imageset
    make_sequence = mocker.patch.object(ExperimentListDict, "_make_sequence")
    # Ensure that if make_sequence is called more than once it returns a different
    # value each time
    make_sequence.side_effect = lambda *args, **kwargs: mocker.MagicMock()

    # Decode the dict to get a new experiment list
    experiments2 = ExperimentListDict(d).decode()

    # This function should only be called once per imageset
    make_sequence.assert_called_once()

    # Verify that this experiment is as we expect
    assert len(experiments2) == 3
    assert len(experiments2.imagesets()) == 1
    assert len(experiments2.goniometers()) == 1
    assert len(experiments2.detectors()) == 1
    assert len(experiments2.beams()) == 1
    assert len(experiments2.scans()) == 3
    for expt in experiments2:
        assert expt.imageset is experiments2.imagesets()[0]


def test_experimentlist_change_basis(dials_data):
    experiments = ExperimentList()
    for i in range(4):
        experiments.extend(
            ExperimentList.from_file(
                dials_data("vmxi_proteinase_k_sweeps", pathlib=True)
                / ("experiments_%i.expt" % i),
                check_format=False,
            )
        )
    reindexed_uc = (68.368, 103.968, 68.368, 90.000, 90.000, 90.000)
    reindexed_sg = sgtbx.space_group_info("P 4 2 2 (b,c,a)").group()
    cb_op = sgtbx.change_of_basis_op("-a,-c,-b")
    for cb_op in (cb_op, [cb_op] * len(experiments)):
        expts_rdx = experiments.change_basis(cb_op)
        for expt in expts_rdx:
            assert expt.crystal.get_unit_cell().parameters() == pytest.approx(
                reindexed_uc, abs=0.1
            )
            assert expt.crystal.get_space_group() == reindexed_sg

    experiments.change_basis(cb_op, in_place=True)
    for expt in experiments:
        assert expt.crystal.get_unit_cell().parameters() == pytest.approx(
            reindexed_uc, abs=0.1
        )
        assert expt.crystal.get_space_group() == reindexed_sg

    with pytest.raises(AssertionError):
        experiments.change_basis([cb_op, cb_op])


def test_path_iterator(monkeypatch):
    """Test the pathname iterator that avoids excessive file calls"""

    @classmethod
    def _fake_open_file(cls, name):
        """Mock replacement for Format's open_file"""
        if name in ("a", "b", os.path.join("dir", "c"), os.path.join("dir", "d"), "e"):
            return mock.Mock()
        elif name.startswith("dir"):
            err = IOError()
            err.errno = errno.EISDIR
            # raise IOError(errno=errno.EISDIR)
            raise err
        assert False

    # Path the lookup of files
    listdir = mock.Mock(return_value=["c", "dir2", "d"])
    monkeypatch.setattr(os, "listdir", listdir)
    # Replace Format.open_file with a tame version
    monkeypatch.setattr(Format, "open_file", _fake_open_file)

    it = dxtbx.model.experiment_list._openingpathiterator(["a", "b", "dir", "e"])
    assert list(it) == [
        "a",
        "b",
        os.path.join("dir", "c"),
        os.path.join("dir", "d"),
        "e",
    ]
    listdir.assert_called_once_with("dir")

    # Test that the list is sorted
    it = dxtbx.model.experiment_list._openingpathiterator(["e", "a", "b", "dir"])
    assert list(it) == [
        "a",
        "b",
        os.path.join("dir", "c"),
        os.path.join("dir", "d"),
        "e",
    ]


def test_extract_metadata_record():
    """Make sure we can read a metadataobject from a format instance"""
    fmt = mock.MagicMock()
    fmt.get_image_file.return_value = "filename_000.cbf"
    fmt.get_scan.return_value = None
    record = dxtbx.model.experiment_list.ImageMetadataRecord.from_format(fmt)
    assert record.beam is fmt.get_beam()
    assert record.detector is fmt.get_detector()
    assert record.goniometer is fmt.get_goniometer()
    assert record.scan is None
    assert record.index is None


def _equal_but_not_same(thing):
    object_1 = (thing,)
    object_2 = (thing,)
    assert object_1 == object_2
    assert object_1 is not object_2
    return object_1, object_2


def test_merge_metadata_record():
    """Test that merging metadata records works correctly"""
    beam_a, beam_b = _equal_but_not_same("beam")
    detector_a, detector_b = _equal_but_not_same("detector")
    gonio_a, gonio_b = _equal_but_not_same("goniometer")

    a = dxtbx.model.experiment_list.ImageMetadataRecord(
        beam=beam_a, detector=detector_a, goniometer=gonio_a
    )
    b = dxtbx.model.experiment_list.ImageMetadataRecord(
        beam=beam_b, detector=detector_b, goniometer=gonio_b
    )
    pre_hash = hash(a)
    assert a.beam is not b.beam
    assert a.detector is not b.detector
    assert a.goniometer is not b.goniometer
    # This should do something
    assert b.merge_metadata_from(a)
    assert hash(a) == pre_hash, "a changed after merge"
    # Make sure metadata was merged
    assert a.beam is b.beam
    assert a.detector is b.detector
    assert a.goniometer is b.goniometer
    # This should NOT do something
    assert not a.merge_metadata_from(a)
    assert hash(a) == pre_hash


def test_merge_all_metadata():
    """Test that merging metadata over a whole list of records works"""
    beam_a, beam_b = _equal_but_not_same("beam")
    gonio_a, gonio_b = _equal_but_not_same("goniometer")

    a = dxtbx.model.experiment_list.ImageMetadataRecord(
        beam=beam_a, detector=object(), goniometer=gonio_a
    )
    b = dxtbx.model.experiment_list.ImageMetadataRecord(
        beam=beam_b, detector=object(), goniometer=gonio_b
    )
    records = [a, b]
    dxtbx.model.experiment_list._merge_model_metadata(records)
    assert a.beam is b.beam
    assert a.goniometer is b.goniometer
    assert a.detector is not b.detector


def test_merge_scan():
    """Test merging logic of scans"""
    scanA = mock.Mock(spec=Scan)
    scanB = mock.Mock(spec=Scan)
    recordA = mock.Mock(scan=scanA)
    recordB = mock.Mock(
        scan=scanB,
        beam=recordA.beam,
        detector=recordA.detector,
        goniometer=recordA.goniometer,
    )
    result = dxtbx.model.experiment_list._merge_scans([recordA, recordB])
    assert result == [recordA]
    scanA.append.assert_called_once_with(scanB)

    # Change some metadata in recordB so it doesn't match
    scanA.reset_mock()
    recordB.beam = mock.Mock()
    assert dxtbx.model.experiment_list._merge_scans([recordA, recordB]) == [
        recordA,
        recordB,
    ]


def test_groupby_template_none():
    Fo = collections.namedtuple("Fo", ["template"])
    objs = [Fo(1), Fo(2), Fo(2), Fo(None), Fo(None), Fo("something")]
    result = list(dxtbx.model.experiment_list._groupby_template_is_none(objs))
    assert result == [
        [Fo(1)],
        [Fo(2)],
        [Fo(2)],
        [Fo(None), Fo(None)],
        [Fo("something")],
    ]


def test_iterate_with_previous():
    sample = list(range(5))
    assert list(dxtbx.model.experiment_list._iterate_with_previous(sample)) == [
        (None, 0),
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
    ]


def test_create_single_sequence(single_sequence_filenames):
    experiments = ExperimentListFactory.from_filenames(single_sequence_filenames)
    assert len(experiments) == 1
    imagesets = experiments.imagesets()
    assert imagesets[0].get_format_class()
    assert len(imagesets) == 1
    assert len(imagesets[0]) == 9


def test_create_multiple_sequences(multiple_sequence_filenames):
    experiments = ExperimentListFactory.from_filenames(multiple_sequence_filenames)
    assert len(experiments) == 2
    imagesets = experiments.imagesets()
    assert len(imagesets) == 2
    assert imagesets[0].get_format_class()
    assert imagesets[1].get_format_class()
    assert len(imagesets[0]) == 3
    assert len(imagesets[1]) == 3


def test_create_multiple_blocks(multiple_block_filenames):
    experiments = ExperimentListFactory.from_filenames(multiple_block_filenames)
    assert len(experiments) == 24
    imagesets = experiments.imagesets()
    assert len(imagesets) == 24
    assert [len(im) for im in imagesets] == [9] + [1] * 23


def test_from_null_sequence():
    filenames = ["template_%2d.cbf" % (i + 1) for i in range(0, 10)]
    sequence = Format.get_imageset(
        filenames,
        beam=Beam((0, 0, 1)),
        detector=Detector(),
        goniometer=Goniometer((1, 0, 0)),
        scan=Scan((1, 10), (0, 0.1)),
    )

    # Create the experiments
    experiments = ExperimentListFactory.from_sequence_and_crystal(
        sequence, crystal=None
    )
    assert len(experiments) == 1
    imagesets = experiments.imagesets()
    assert imagesets[0].get_format_class()
    assert len(imagesets) == 1
    assert imagesets[0].get_beam() == sequence.get_beam()
    assert imagesets[0].get_detector() == sequence.get_detector()
    assert imagesets[0].get_goniometer() == sequence.get_goniometer()
    assert imagesets[0].get_scan() == sequence.get_scan()


def test_from_templates(dials_data):
    template = dials_data("insulin", pathlib=True) / "insulin_1_###.img"
    expts = ExperimentList.from_templates([template])
    assert len(expts) == 1
    assert expts[0].imageset.get_template() == str(template)
    assert len(expts[0].imageset) == 45


def test_experiment_list_all():
    experiments = ExperimentList()
    for i in range(3):
        experiments.append(Experiment())

    assert experiments.all_stills()
    experiments[0].goniometer = Goniometer()
    assert experiments.all_stills()
    experiments[1].goniometer = Goniometer()
    experiments[2].goniometer = Goniometer()
    assert experiments.all_stills()

    experiments[0].beam = BeamFactory.make_polychromatic_beam(
        direction=(0, 0, -1),
        sample_to_source_distance=(100),
        probe=Probe.xray,
        wavelength_range=(1, 10),
    )
    assert not experiments.all_stills()
    experiments[1].beam = BeamFactory.make_polychromatic_beam(
        direction=(0, 0, -1),
        sample_to_source_distance=(100),
        probe=Probe.xray,
        wavelength_range=(1, 10),
    )
    experiments[2].beam = BeamFactory.make_polychromatic_beam(
        direction=(0, 0, -1),
        sample_to_source_distance=(100),
        probe=Probe.xray,
        wavelength_range=(1, 10),
    )
    assert experiments.all_laue()

    experiments[0].beam = Beam()
    assert not experiments.all_laue()
    experiments[1].beam = Beam()
    experiments[2].beam = Beam()
    assert experiments.all_stills()

    experiments[0].scan = Scan((1, 1000), (0, 0.05))
    assert not experiments.all_stills()
    experiments[1].scan = Scan((1, 1000), (0, 0.05))
    experiments[2].scan = Scan((1, 1000), (0, 0.05))
    assert experiments.all_rotations()

    experiments[0].scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"time_of_flight": list(range(10))}
    )
    assert not experiments.all_rotations()
    experiments[1].scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"time_of_flight": list(range(10))}
    )
    experiments[2].scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"time_of_flight": list(range(10))}
    )
    assert experiments.all_tof()

    experiments[0].scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"other_property": list(range(10))}
    )
    assert not experiments.all_tof()
    experiments[1].scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"other_property": list(range(10))}
    )
    experiments[2].scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"other_property": list(range(10))}
    )
    assert experiments.all_stills()
