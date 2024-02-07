from __future__ import annotations

import pickle

import pytest

import cctbx.array_family.flex as flex

from dxtbx.model import Scan, ScanFactory


def test_scan_to_string_does_not_crash_on_empty_scan():
    print(Scan())


def test_scan_wrap_around_zero():
    filenames = [f"foobar_{n:03d}.cbf" for n in range(350, 370)]
    starts = [r % 360 for r in range(350, 370)]
    ends = [(r + 1) % 360 for r in range(350, 370)]
    scans = [
        ScanFactory.single_file(f, [0.1], s, e - s, 0)
        for f, s, e in zip(filenames, starts, ends)
    ]
    s0 = scans[0]
    for s in scans[1:]:
        assert s.get_oscillation()[1] == 1
        s0 += s
    assert s0.get_oscillation() == pytest.approx((350, 1))
    assert s0.get_image_range() == (350, 369)


def test_make_scan_from_properties():
    image_range = (1, 10)

    properties = {}
    scan = ScanFactory.make_scan_from_properties(
        image_range=image_range, properties=properties
    )
    assert scan.get_properties() == {}

    properties = {
        "test_int": tuple(range(10)),
        "test_float": tuple([float(i * 2) for i in range(10)]),
        "test_bool": tuple([True for i in range(10)]),
        "test_string": tuple([f"test_{i}" for i in range(10)]),
        "test_vec3_double": tuple([(1.0, 1.0, 1.0) for i in range(10)]),
        "test_vec2_double": tuple([(2.0, 2.0) for i in range(10)]),
    }

    scan = ScanFactory.make_scan_from_properties(
        image_range=image_range, properties=properties
    )

    result = scan.get_properties()
    for key, value in properties.items():
        assert key in result
        assert len(value) == len(result[key])
        for i in range(len(value)):
            assert result[key][i] == pytest.approx(value[i])


def test_set_get_properties():

    image_range = (1, 10)
    scan = ScanFactory.make_scan_from_properties(image_range=image_range, properties={})

    test_int = flex.int(10, 1)
    test_float = flex.double(10, 1.0)
    test_bool = flex.bool(10, True)
    test_string = flex.std_string(10, "Test")
    test_vec3_double = flex.vec3_double(10, (1.0, 1.0, 1.0))
    test_vec2_double = flex.vec2_double(10, (2.0, 2.0))

    scan.set_property("test_int", test_int)
    result = scan.get_property("test_int")
    for idx, i in enumerate(test_int):
        assert result[idx] == pytest.approx(i)

    scan.set_property("test_float", test_float)
    result = scan.get_property("test_float")
    for idx, i in enumerate(test_float):
        assert result[idx] == pytest.approx(i)

    scan.set_property("test_bool", test_bool)
    result = scan.get_property("test_bool")
    for idx, i in enumerate(test_bool):
        assert result[idx] == pytest.approx(i)

    scan.set_property("test_string", test_string)
    result = scan.get_property("test_string")
    for idx, i in enumerate(test_string):
        assert result[idx] == pytest.approx(i)

    scan.set_property("test_vec3_double", test_vec3_double)
    result = scan.get_property("test_vec3_double")
    for idx, i in enumerate(test_vec3_double):
        assert result[idx] == pytest.approx(i)

    scan.set_property("test_vec2_double", test_vec2_double)
    result = scan.get_property("test_vec2_double")
    for idx, i in enumerate(test_vec2_double):
        assert result[idx] == pytest.approx(i)

    # Test incompatible size with property table

    test_int = flex.int(5, 1)
    with pytest.raises(RuntimeError):
        scan.set_property("test_int", test_int)

    test_float = flex.double(5, 1.0)
    with pytest.raises(RuntimeError):
        scan.set_property("test_float", test_float)

    test_bool = flex.bool(5, True)
    with pytest.raises(RuntimeError):
        scan.set_property("test_bool", test_bool)

    test_string = flex.std_string(5, "Test")
    with pytest.raises(RuntimeError):
        scan.set_property("test_string", test_string)

    test_vec3_double = flex.vec3_double(5, (1.0, 1.0, 1.0))
    with pytest.raises(RuntimeError):
        scan.set_property("test_vec3_double", test_vec3_double)

    test_vec2_double = flex.vec2_double(5, (2.0, 2.0))
    with pytest.raises(RuntimeError):
        scan.set_property("test_vec2_double", test_vec2_double)

    test_int = flex.int(15, 1)
    with pytest.raises(RuntimeError):
        scan.set_property("test_int", test_int)

    test_float = flex.double(15, 1.0)
    with pytest.raises(RuntimeError):
        scan.set_property("test_float", test_float)

    test_bool = flex.bool(15, True)
    with pytest.raises(RuntimeError):
        scan.set_property("test_bool", test_bool)

    test_string = flex.std_string(15, "Test")
    with pytest.raises(RuntimeError):
        scan.set_property("test_string", test_string)

    test_vec3_double = flex.vec3_double(15, (1.0, 1.0, 1.0))
    with pytest.raises(RuntimeError):
        scan.set_property("test_vec3_double", test_vec3_double)

    test_vec2_double = flex.vec2_double(15, (2.0, 2.0))
    with pytest.raises(RuntimeError):
        scan.set_property("test_vec2_double", test_vec2_double)


def test_scan_properties_pickle():
    image_range = (1, 10)
    properties = {
        "test_int": tuple(range(10)),
        "test_float": tuple([float(i * 2) for i in range(10)]),
        "test_bool": tuple([True for i in range(10)]),
        "test_string": tuple([f"test_{i}" for i in range(10)]),
        "test_vec3_double": tuple([(1.0, 1.0, 1.0) for i in range(10)]),
        "test_vec2_double": tuple([(2.0, 2.0) for i in range(10)]),
    }

    scan = ScanFactory.make_scan_from_properties(
        image_range=image_range, properties=properties
    )

    obj = pickle.dumps(scan)
    scan2 = pickle.loads(obj)
    assert scan == scan2

    result = scan.get_properties()
    for key, value in properties.items():
        assert key in result
        assert len(value) == len(result[key])
        for i in range(len(value)):
            assert result[key][i] == pytest.approx(value[i])


def test_scan_properties_to_dict():
    image_range = (1, 10)
    properties = {
        "test_int": tuple(range(10)),
        "test_float": tuple([float(i * 2) for i in range(10)]),
        "test_bool": tuple([True for i in range(10)]),
        "test_string": tuple([f"test_{i}" for i in range(10)]),
        "test_vec3_double": tuple([(1.0, 1.0, 1.0) for i in range(10)]),
        "test_vec2_double": tuple([(2.0, 2.0) for i in range(10)]),
    }

    scan = ScanFactory.make_scan_from_properties(
        image_range=image_range, properties=properties
    )

    expected_result = {
        "test_bool": (True, True, True, True, True, True, True, True, True, True),
        "test_float": (0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0),
        "test_int": (0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
        "test_string": (
            "test_0",
            "test_1",
            "test_2",
            "test_3",
            "test_4",
            "test_5",
            "test_6",
            "test_7",
            "test_8",
            "test_9",
        ),
        "test_vec2_double": (
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
        ),
        "test_vec3_double": (
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
        ),
    }

    result = scan.to_dict()["properties"]
    for key, value in expected_result.items():
        assert key in result
        assert len(value) == len(result[key])
        for i in range(len(value)):
            assert result[key][i] == pytest.approx(value[i])


def test_scan_properties_equivalence():

    int_diff = 1
    double_diff = 1e-6

    image_range = (1, 10)
    s1 = ScanFactory.make_scan_from_properties(image_range=image_range, properties={})
    s2 = ScanFactory.make_scan_from_properties(image_range=image_range, properties={})

    # Empty table

    assert s1 == s2

    # int columns

    test_int = flex.int(10, 1)
    test_int2 = flex.int(10, 1)
    test_int2[0] += int_diff

    s1.set_property("test_int", test_int)

    assert s1 != s2

    s2.set_property("test_int", test_int2)

    assert s1 != s2

    test_int2[0] -= int_diff

    s2.set_property("test_int", test_int2)

    assert s1 == s2

    # double columns

    test_double = flex.double(10, 1.0)
    test_double2 = flex.double(10, 1.0)
    test_double2[0] += double_diff

    s1.set_property("test_double", test_double)

    assert s1 != s2

    s2.set_property("test_double", test_double2)

    assert s1 != s2

    test_double2[0] -= double_diff

    s2.set_property("test_double", test_double2)

    assert s1 == s2

    # bool columns

    test_bool = flex.bool(10, True)
    test_bool2 = flex.bool(10, True)
    test_bool2[0] = False

    s1.set_property("test_bool", test_bool)

    assert s1 != s2

    s2.set_property("test_bool", test_bool2)

    assert s1 != s2

    test_bool2[0] = True

    s2.set_property("test_bool", test_bool2)

    assert s1 == s2

    # string columns

    test_string = flex.std_string(10, "Test")
    test_string2 = flex.std_string(10, "Test")
    test_string2[0] = "Test2"

    s1.set_property("test_string", test_string)

    assert s1 != s2

    s2.set_property("test_string", test_string2)

    assert s1 != s2

    test_string2[0] = "Test"

    s2.set_property("test_string", test_string2)

    assert s1 == s2

    #  vec3<double> columns

    test_vec3_double = flex.vec3_double(10, (1.0, 1.0, 1.0))
    test_vec3_double2 = flex.vec3_double(10, (1.0, 1.0, 1.0))
    test_vec3_double2[0] = (1.0 + double_diff, 1.0 + double_diff, 1.0 + double_diff)

    s1.set_property("test_vec3_double", test_vec3_double)

    assert s1 != s2

    s2.set_property("test_vec3_double", test_vec3_double2)

    assert s1 != s2
    vec = test_vec3_double2[0]
    test_vec3_double2[0] = (
        vec[0] - double_diff,
        vec[1] - double_diff,
        vec[2] - double_diff,
    )

    s2.set_property("test_vec3_double", test_vec3_double2)

    assert s1 == s2

    #  vec2<double> columns

    test_vec2_double = flex.vec2_double(10, (2.0, 2.0))
    test_vec2_double2 = flex.vec2_double(10, (2.0, 2.0))
    test_vec2_double2[0] = (2.0 + double_diff, 2.0 + double_diff)

    s1.set_property("test_vec2_double", test_vec2_double)

    assert s1 != s2

    s2.set_property("test_vec2_double", test_vec2_double2)

    assert s1 != s2

    vec = test_vec2_double2[0]
    test_vec2_double2[0] = (vec[0] - double_diff, vec[1] - double_diff)

    s2.set_property("test_vec2_double", test_vec2_double2)

    assert s1 == s2


def test_scan_constant_oscillation_width():

    image_range = (1, 10)
    scan = ScanFactory.make_scan_from_properties(image_range=image_range, properties={})

    oscillation = flex.double(10, 1.0)
    scan.set_property("oscillation", oscillation)
    oscillation[3] = 2.0
    with pytest.raises(RuntimeError):
        scan.set_property("oscillation", oscillation)

    oscillation = flex.double(10, 1.0)
    properties = {"oscillation": oscillation}
    scan.set_properties(properties)
    properties["oscillation"][3] = 2.0
    with pytest.raises(RuntimeError):
        scan.set_properties(properties)

    oscillation = flex.double(10, 1.0)
    properties = {"oscillation": oscillation}
    scan = ScanFactory.make_scan_from_properties(
        image_range=image_range, properties=properties
    )
    properties["oscillation"][3] = 2.0
    with pytest.raises(RuntimeError):
        scan = ScanFactory.make_scan_from_properties(
            image_range=image_range, properties=properties
        )


def test_print_scan():

    image_range = (1, 10)
    scan = ScanFactory.make_scan_from_properties(image_range=image_range, properties={})

    expected_scan_string = (
        "Scan:\n    number of images:   10\n    image range:   {1,10}\n"
    )
    assert scan.__str__() == expected_scan_string

    # Properties with unique output

    scan.set_oscillation((0.0, 0.5))
    expected_scan_string = "Scan:\n    number of images:   10\n    image range:   {1,10}\n    oscillation:   {0,0.5}\n"
    assert scan.__str__() == expected_scan_string

    scan.set_property("epochs", flex.double(10))
    expected_scan_string = "Scan:\n    number of images:   10\n    image range:   {1,10}\n    epoch:    0\n    oscillation:   {0,0.5}\n"
    assert scan.__str__() == expected_scan_string

    scan.set_property("exposure_time", flex.double(10))
    expected_scan_string = "Scan:\n    number of images:   10\n    image range:   {1,10}\n    epoch:    0\n    exposure time:    0\n    oscillation:   {0,0.5}\n"
    assert scan.__str__() == expected_scan_string

    # Generic properties

    properties = {
        "test_int": tuple(range(10)),
        "test_float": tuple([float(i * 2) for i in range(10)]),
        "test_bool": tuple([True for i in range(10)]),
        "test_string": tuple([f"test_{i}" for i in range(10)]),
        "test_vec3_double": tuple([(1.0, 1.0, 1.0) for i in range(10)]),
        "test_vec2_double": tuple([(2.0, 2.0) for i in range(10)]),
    }

    scan = ScanFactory.make_scan_from_properties(
        image_range=image_range, properties=properties
    )
    expected_scan_string = "Scan:\n    number of images:   10\n    image range:   {1,10}\n    test_bool:    1 - 1\n    test_float:    0 - 18\n    test_int:    0 - 9\n    test_string:    test_0 - test_9\n    test_vec2_double:    {2,2} - {2,2}\n    test_vec3_double:    {1,1,1} - {1,1,1}\n"
    assert scan.__str__() == expected_scan_string


def test_scan_is_still():
    scan = Scan()
    assert scan.is_still()
    scan = Scan((1, 1000), (0, 0.05))
    assert not scan.is_still()
    # Specifically test the bug from dxtbx#4 triggered by ending on 0°
    scan = Scan((1, 1800), (-90, 0.05))
    assert not scan.is_still()
    scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"time_of_flight": list(range(10))}
    )
    assert not scan.is_still()
    scan = ScanFactory.make_scan_from_properties(
        (1, 10), properties={"other_property": list(range(10))}
    )
    assert scan.is_still()
