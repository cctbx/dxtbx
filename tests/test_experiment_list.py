from dxtbx.model.experiment_list import ExperimentListFactory


def test_experiment_list_from_filenames(dials_data):
    filenames = [
        dials_data("centroid_test_data").join("centroid_%04d.cbf" % i).strpath
        for i in range(1, 10)
    ]
    experiments = ExperimentListFactory.from_filenames(filenames)

    assert len(experiments) == 1
    expt = experiments[0]
    assert len(expt.imageset) == 9
    assert expt.index == 0


def test_experiment_list_from_filenames_missing_one(dials_data):
    filenames = [
        dials_data("centroid_test_data").join("centroid_%04d.cbf" % i).strpath
        for i in range(1, 10)
    ]
    filenames.pop(4)
    experiments = ExperimentListFactory.from_filenames(filenames)

    assert len(experiments) == 2
    for j, expt in enumerate(experiments):
        assert len(expt.imageset) == 4
        assert expt.index == j
