from __future__ import annotations

from os.path import join

import pytest

from dxtbx.format.FormatISISSXD import FormatISISSXD


@pytest.fixture(scope="session")
def nacl(dials_data):
    location = dials_data("isis_sxd_example_data", pathlib=True)
    nacl_filename = join(location, "sxd_nacl_run.nxs")
    return nacl_filename


def test_import(nacl):
    assert FormatISISSXD.understand(nacl)
    fmt = FormatISISSXD(nacl)

    assert fmt.get_experiment_title() == "NaCl sphere 6mm diameter RT j:14,14"
    assert fmt.get_experiment_run_number() == 33298
    assert (
        fmt.get_experiment_description()
        == "NaCl sphere 6mm diameter RT j:14,14 (33298)"
    )

    goniometer = fmt.get_goniometer()
    assert goniometer.to_dict() == {
        "rotation_axis": (0.0, 1.0, 0.0),
        "fixed_rotation": (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        "setting_rotation": (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    }

    beam = fmt.get_beam()
    assert beam.to_dict() == {
        "__id__": "polychromatic",
        "direction": (0.0, 0.0, -1.0),
        "divergence": 0.0,
        "sigma_divergence": 0.0,
        "polarization_normal": (0.0, 1.0, 0.0),
        "polarization_fraction": 0.5,
        "flux": 0.0,
        "transmission": 1.0,
        "probe": "neutron",
        "sample_to_source_distance": 8300.0,
        "wavelength_range": (0.2, 10.0),
    }

    detector = fmt.get_detector()
    det_dict = detector.to_dict()
    for p in det_dict["panels"]:
        p["fast_axis"] = tuple(round(elt, 10) for elt in p["fast_axis"])
        p["slow_axis"] = tuple(round(elt, 10) for elt in p["slow_axis"])
    assert det_dict == {
        "panels": [
            {
                "name": "01",
                "type": "SENSOR_PAD",
                "fast_axis": (0.0, -1.0, 0.0),
                "slow_axis": (0.7931070767, 0.0, 0.6090822317),
                "origin": (60.81, 96.0, -236.946),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (-64, -128)),
            },
            {
                "name": "02",
                "type": "SENSOR_PAD",
                "fast_axis": (-0.0, -0.0, -1.0),
                "slow_axis": (0.0, 1.0, 0.0),
                "origin": (224.999, -96.0, 96.0),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (0, -128)),
            },
            {
                "name": "03",
                "type": "SENSOR_PAD",
                "fast_axis": (0.7931070767, -0.0, -0.6090822317),
                "slow_axis": (0.0, 1.0, 0.0),
                "origin": (60.809, -96.0, 236.945),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (64, -128)),
            },
            {
                "name": "04",
                "type": "SENSOR_PAD",
                "fast_axis": (0.7878424473, -0.0, 0.6158768369),
                "slow_axis": (0.0, 1.0, 0.0),
                "origin": (-214.172, -96.0, 118.198),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (64, 128)),
            },
            {
                "name": "05",
                "type": "SENSOR_PAD",
                "fast_axis": (-0.0, -0.0, 1.0),
                "slow_axis": (0.0, 1.0, 0.0),
                "origin": (-224.999, -96.0, -96.0),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (0, 128)),
            },
            {
                "name": "06",
                "type": "SENSOR_PAD",
                "fast_axis": (-0.7931070767, -0.0, 0.6090822317),
                "slow_axis": (0.0, 1.0, 0.0),
                "origin": (-60.809, -96.0, -236.945),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (-64, 128)),
            },
            {
                "name": "07",
                "type": "SENSOR_PAD",
                "fast_axis": (0.0, -0.0, -1.0),
                "slow_axis": (0.6950048651, 0.7190050331, -0.0),
                "origin": (127.534, -256.614, 96.0),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (0, -64)),
            },
            {
                "name": "08",
                "type": "SENSOR_PAD",
                "fast_axis": (1.0, -0.0, -0.0),
                "slow_axis": (0.0, 0.7190050331, 0.6950048651),
                "origin": (-96.0, -256.614, 127.534),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (64, 0)),
            },
            {
                "name": "09",
                "type": "SENSOR_PAD",
                "fast_axis": (-0.0, -0.0, 1.0),
                "slow_axis": (-0.7071067812, 0.7071067812, -0.0),
                "origin": (-123.036, -258.801, -96.0),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (0, 64)),
            },
            {
                "name": "10",
                "type": "SENSOR_PAD",
                "fast_axis": (-1.0, -0.0, -0.0),
                "slow_axis": (0.0, 0.7190050331, -0.6950048651),
                "origin": (96.0, -256.614, -127.534),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (-64, 0)),
            },
            {
                "name": "11",
                "type": "SENSOR_PAD",
                "fast_axis": (-1.0, -0.0, -0.0),
                "slow_axis": (-0.0, 0.0, -1.0),
                "origin": (96.0, -278.0, 96.0),
                "raw_image_offset": (0, 0),
                "image_size": (64, 64),
                "pixel_size": (3.0, 3.0),
                "trusted_range": (-1.0, 100000.0),
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                "projection_2d": ((1, 0, 0, 1), (0, 0)),
            },
        ],
        "hierarchy": {
            "name": "",
            "type": "",
            "fast_axis": (1.0, 0.0, 0.0),
            "slow_axis": (0.0, 1.0, 0.0),
            "origin": (0.0, 0.0, 0.0),
            "raw_image_offset": (0, 0),
            "image_size": (0, 0),
            "pixel_size": (0.0, 0.0),
            "trusted_range": (0.0, 0.0),
            "thickness": 0.0,
            "material": "",
            "mu": 0.0,
            "identifier": "",
            "mask": [],
            "gain": 1.0,
            "pedestal": 0.0,
            "px_mm_strategy": {"type": "SimplePxMmStrategy"},
            "children": [
                {"panel": 0},
                {"panel": 1},
                {"panel": 2},
                {"panel": 3},
                {"panel": 4},
                {"panel": 5},
                {"panel": 6},
                {"panel": 7},
                {"panel": 8},
                {"panel": 9},
                {"panel": 10},
            ],
        },
    }

    scan = fmt.get_scan()
    assert len(scan.get_properties())
    assert scan.has_property("time_of_flight")
    assert scan.get_image_range() == (1, 1821)
    assert scan.get_property("time_of_flight")[0] == pytest.approx(500.5)
