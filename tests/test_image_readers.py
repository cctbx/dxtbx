from __future__ import absolute_import, division, print_function

import binascii
import os
from builtins import range

import numpy
import pytest

from boost.python import streambuf
from cbflib_adaptbx import uncompress
from scitbx.array_family import flex

import dxtbx.ext
import dxtbx.format.FormatCBF
import dxtbx.tests.imagelist
import pycbf
from dxtbx.format.FormatSMV import FormatSMV
from dxtbx.format.FormatTIFFHelpers import read_basic_tiff_header
from dxtbx.format.image import (
    CBFReader,
    HDF5Reader,
    SMVReader,
    TIFFReader,
    cbf_read_buffer,
)
from dxtbx.model.detector import DetectorFactory


def read_smv_image(image_file):
    header_size, header_dictionary = FormatSMV.get_smv_header(image_file)

    with open(image_file, "rb") as f:
        f.seek(header_size)

        big_endian = header_dictionary["BYTE_ORDER"] == "big_endian"

        image_size = (int(header_dictionary["SIZE1"]), int(header_dictionary["SIZE2"]))

        if big_endian == dxtbx.ext.is_big_endian():
            raw_data = dxtbx.ext.read_uint16(
                streambuf(f), int(image_size[0] * image_size[1])
            )
        else:
            raw_data = dxtbx.ext.read_uint16_bs(
                streambuf(f), int(image_size[0] * image_size[1])
            )

    raw_data.reshape(flex.grid(image_size[1], image_size[0]))

    return raw_data


def get_tiff_header(image_file):
    width, height, depth, header, order = read_basic_tiff_header(image_file)

    with open(image_file, "rb") as fh:
        header_bytes = fh.read(header)

    return width, height, depth // 8, order, header_bytes


def read_tiff_image(image_file):
    # currently have no non-little-endian machines...
    width, height, depth, order, header_bytes = get_tiff_header(image_file)
    image_size = (width, height)
    header_size = 4096

    with open(image_file, "rb") as fh:
        fh.seek(header_size)
        raw_data = dxtbx.ext.read_uint16(
            streambuf(fh), int(image_size[0] * image_size[1])
        )
    raw_data.reshape(flex.grid(image_size[1], image_size[0]))

    return raw_data


def read_cbf_image(cbf_image):
    start_tag = binascii.unhexlify("0c1a04d5")

    with open(cbf_image, "rb") as fh:
        data = fh.read()
    data_offset = data.find(start_tag) + 4
    cbf_header = dxtbx.format.FormatCBF.FormatCBF._parse_cbf_header(
        data[: data_offset - 4].decode("ascii", "ignore")
    )

    pixel_values = uncompress(
        packed=data[data_offset : data_offset + cbf_header["size"]],
        fast=cbf_header["fast"],
        slow=cbf_header["slow"],
    )

    return pixel_values


def read_multitile_cbf_image(cbf_image):
    raw_data = []
    cbf = pycbf.cbf_handle_struct()
    cbf.read_widefile(cbf_image.encode(), pycbf.MSG_DIGEST)
    cbf.find_category(b"array_structure")
    cbf.find_column(b"encoding_type")
    cbf.select_row(0)
    types = []
    for i in range(cbf.count_rows()):
        types.append(cbf.get_value())
        cbf.next_row()
    assert len(types) == cbf.count_rows()

    # read the data
    data = {}
    cbf.find_category(b"array_data")
    for i in range(cbf.count_rows()):
        cbf.find_column(b"array_id")
        name = cbf.get_value()

        cbf.find_column(b"data")
        assert cbf.get_typeofvalue().find(b"bnry") > -1

        if types[i] == b"signed 32-bit integer":
            array_string = cbf.get_integerarray_as_string()
            array = flex.int(numpy.frombuffer(array_string, numpy.int32))
            parameters = cbf.get_integerarrayparameters_wdims_fs()
            array_size = (parameters[11], parameters[10], parameters[9])
        elif types[i] == b"signed 64-bit real IEEE":
            array_string = cbf.get_realarray_as_string()
            array = flex.double(numpy.frombuffer(array_string, numpy.float))
            parameters = cbf.get_realarrayparameters_wdims_fs()
            array_size = (parameters[7], parameters[6], parameters[5])
        else:
            return None  # type not supported

        array.reshape(flex.grid(*array_size))
        data[name] = array
        cbf.next_row()

    # extract the data for each panel
    try:
        cbf.find_category(b"array_structure_list_section")
        has_sections = True
    except Exception:
        has_sections = False
    if has_sections:
        section_shapes = {}
        for i in range(cbf.count_rows()):
            cbf.find_column(b"id")
            section_name = cbf.get_value()
            if section_name not in section_shapes:
                section_shapes[section_name] = {}
            cbf.find_column(b"array_id")
            if "array_id" in section_shapes[section_name]:
                assert section_shapes[section_name]["array_id"] == cbf.get_value()
            else:
                section_shapes[section_name]["array_id"] = cbf.get_value()
            cbf.find_column(b"index")
            axis_index = int(cbf.get_value()) - 1
            cbf.find_column(b"start")
            axis_start = int(cbf.get_value()) - 1
            cbf.find_column(b"end")
            axis_end = int(cbf.get_value())

            section_shapes[section_name][axis_index] = slice(axis_start, axis_end)
            cbf.next_row()

        for section_name in sorted(section_shapes):
            section_shape = section_shapes[section_name]
            section = data[section_shape["array_id"]][
                section_shape[2], section_shape[1], section_shape[0]
            ]
            section.reshape(flex.grid(section.focus()[-2], section.focus()[-1]))
            raw_data.append(section)
    else:
        for key in sorted(data):
            data[key].reshape(flex.grid(data[key].focus()[-2], data[key].focus()[-1]))
            raw_data.append(data[key])

    return tuple(raw_data)


@pytest.mark.parametrize(
    "smv_image",
    dxtbx.tests.imagelist.smv_images,
    ids=dxtbx.tests.imagelist.smv_image_ids,
)
def test_smv(dials_regression, smv_image):
    filename = os.path.join(dials_regression, smv_image)

    image = SMVReader(filename).image()
    if image.is_double():
        image = image.as_double()
    else:
        image = image.as_int()
    assert image.n_tiles() == 1
    data1 = image.tile(0).data()

    data2 = read_smv_image(filename)

    diff = flex.abs(data1 - data2)
    assert flex.max(diff) < 1e-7


@pytest.mark.parametrize(
    "tiff_image",
    dxtbx.tests.imagelist.tiff_images,
    ids=dxtbx.tests.imagelist.tiff_image_ids,
)
def test_tiff(dials_regression, tiff_image):
    filename = os.path.join(dials_regression, tiff_image)

    image = TIFFReader(filename).image()
    if image.is_double():
        image = image.as_double()
    else:
        image = image.as_int()
    assert image.n_tiles() == 1
    data1 = image.tile(0).data()

    data2 = read_tiff_image(filename)

    diff = flex.abs(data1 - data2)
    assert flex.max(diff) < 1e-7


@pytest.mark.parametrize(
    "cbf_image",
    dxtbx.tests.imagelist.cbf_images,
    ids=dxtbx.tests.imagelist.cbf_image_ids,
)
def test_cbf(dials_regression, cbf_image):
    filename = os.path.join(dials_regression, cbf_image)

    image = CBFReader(filename).image().as_int()
    assert image.n_tiles() == 1

    data1 = image.tile(0).data()
    data2 = read_cbf_image(filename)
    diff = flex.abs(data1 - data2)
    assert flex.max(diff) < 1e-7


@pytest.mark.parametrize(
    "cbf_image",
    dxtbx.tests.imagelist.cbf_multitile_images,
    ids=dxtbx.tests.imagelist.cbf_multitile_image_ids,
)
def test_multitile_cbf(dials_regression, cbf_image):
    filename = os.path.join(dials_regression, cbf_image)

    image = CBFReader(filename).image().as_int()
    data1 = tuple(image.tile(i).data() for i in range(image.n_tiles()))
    data2 = read_multitile_cbf_image(filename)

    assert len(data1) == len(data2)
    for d1, d2 in zip(data1, data2):
        diff = flex.abs(d1 - d2)
        assert flex.max(diff) < 1e-7


@pytest.mark.parametrize(
    "hdf5_image",
    dxtbx.tests.imagelist.hdf5_images,
    ids=dxtbx.tests.imagelist.hdf5_image_ids,
)
def test_hdf5(dials_regression, hdf5_image):
    h5py = pytest.importorskip("h5py")
    # Import after we know h5py is present
    from dxtbx.format.nexus import dataset_as_flex_int

    filename = os.path.join(dials_regression, hdf5_image)
    handle = h5py.File(filename, "r")

    reader = HDF5Reader(handle.id.id, flex.std_string(["/entry/data/data"]))

    image = reader.image(0)

    data1 = image.as_int().tile(0).data()

    dataset = handle["/entry/data/data"]
    N, height, width = dataset.shape
    data2 = dataset_as_flex_int(
        dataset.id.id, (slice(0, 1, 1), slice(0, height, 1), slice(0, width, 1))
    )
    data2.reshape(flex.grid(data2.all()[1:]))

    assert N == len(reader)
    assert data1.all()[0] == data2.all()[0]
    assert data1.all()[1] == data2.all()[1]
    diff = flex.abs(data1 - data2)
    assert flex.max(diff) < 1e-7


def test_cbf_buffer(dials_regression):
    filename = os.path.join(
        dials_regression, "image_examples", "dials-190", "whatev1_01_00001.cbf"
    )
    with open(filename, "rb") as f:
        contents = f.read()

    handle = pycbf.cbf_handle_struct()
    cbf_read_buffer(handle, contents, pycbf.MSG_DIGEST)
    det = DetectorFactory.imgCIF_H(handle, "unknown")
    assert det
