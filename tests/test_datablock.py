from __future__ import annotations

import os
import pickle

import pytest


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


def pickle_then_unpickle(obj):
    return pickle.loads(pickle.dumps(obj))
