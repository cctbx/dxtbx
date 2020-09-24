/*
 * nexus_ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/slice.hpp>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/error.h>
#include <vector>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "bitshuffle.h"

namespace dxtbx { namespace format { namespace boost_python {

  using namespace boost::python;
  template <typename T>
  inline herr_t custom_read(hid_t,
                            hid_t,
                            hid_t,
                            scitbx::af::versa<T, scitbx::af::flex_grid<> >);

  /**
   * Custom read function for integers
   */
  template <>
  inline herr_t custom_read<int>(
    hid_t dataset_id,
    hid_t mem_space_id,
    hid_t file_space_id,
    scitbx::af::versa<int, scitbx::af::flex_grid<> > data) {
    return H5Dread(
      dataset_id, H5T_NATIVE_INT, mem_space_id, file_space_id, H5P_DEFAULT, &data[0]);
  }

  /**
   * Custom read function for floats
   */
  template <>
  inline herr_t custom_read<float>(
    hid_t dataset_id,
    hid_t mem_space_id,
    hid_t file_space_id,
    scitbx::af::versa<float, scitbx::af::flex_grid<> > data) {
    return H5Dread(
      dataset_id, H5T_NATIVE_FLOAT, mem_space_id, file_space_id, H5P_DEFAULT, &data[0]);
  }

  /**
   * Custom read function for doubles
   */
  template <>
  inline herr_t custom_read<double>(
    hid_t dataset_id,
    hid_t mem_space_id,
    hid_t file_space_id,
    scitbx::af::versa<double, scitbx::af::flex_grid<> > data) {
    return H5Dread(dataset_id,
                   H5T_NATIVE_DOUBLE,
                   mem_space_id,
                   file_space_id,
                   H5P_DEFAULT,
                   &data[0]);
  }


  template <typename T>
  inline scitbx::af::versa<T, scitbx::af::flex_grid<> > image_as_flex(
    hid_t dataset_id, const hsize_t index) {
    // index is the offset

    hid_t datatype = H5Dget_type(dataset_id);
    hsize_t datasize = H5Tget_size(datatype);

    // Do a check on the data
    hid_t file_space_id = H5Dget_space(dataset_id);
    std::size_t rank = H5Sget_simple_extent_ndims(file_space_id);
    DXTBX_ASSERT(rank == 3);

    // Now get the dimensions of the data i.e. image size.
    std::vector<hsize_t> dataset_dims(rank);
    H5Sget_simple_extent_dims(file_space_id, &dataset_dims[0], NULL);

    // Create the data array
    // first make the grid
    scitbx::af::flex_grid<>::index_type dims(2);
    dims[0] = dataset_dims[1];
    dims[1] = dataset_dims[2];
    scitbx::af::flex_grid<> grid(dims);
    // now make the array
    scitbx::af::versa<T, scitbx::af::flex_grid<> > data(
      grid, scitbx::af::init_functor_null<T>());

    // get storage size and allocate to a buffer
    hsize_t size;
    hsize_t offset[3];
    offset[0] = index;
    offset[1] = 0;
    offset[2] = 0;

    herr_t status0 = H5Dget_chunk_storage_size(dataset_id, offset, &size);
    DXTBX_ASSERT(status0 >= 0);
    std::cout << "Got chunk storage size" << std::endl;
    char * chunk = (char *)malloc(size);

    // read the raw data to the buffer (chunk)
    uint32_t filter = 0;
    herr_t status1 = H5DOread_chunk(dataset_id, H5P_DEFAULT, offset, &filter, chunk);
    DXTBX_ASSERT(status1 >= 0);
    std::cout << "Read chunk " << std::endl;

    // decompress using bshuf to the data array
    herr_t status2 = bshuf_decompress_lz4(chunk + 12, &data, size, datasize, 0);
    DXTBX_ASSERT(status2 >= 0);
    std::cout << "Decompressed " << std::endl;

    free(chunk);

    return data;
  }


  /**
   * A function to extract data from hdf5 dataset directly into a flex array
   */
  template <typename T>
  inline scitbx::af::versa<T, scitbx::af::flex_grid<> > dataset_as_flex(
    hid_t dataset_id,
    boost::python::tuple selection) {
    // The number of dimensions
    std::size_t ndims = boost::python::len(selection);

    // Get the file space
    hid_t file_space_id = H5Dget_space(dataset_id);
    std::size_t rank = H5Sget_simple_extent_ndims(file_space_id);
    DXTBX_ASSERT(rank == ndims);
    std::vector<hsize_t> dataset_dims(rank);
    H5Sget_simple_extent_dims(file_space_id, &dataset_dims[0], NULL);

    // Create the grid
    scitbx::af::flex_grid<>::index_type dims(ndims);
    std::vector<hsize_t> start(ndims);
    std::vector<hsize_t> count(ndims);
    for (std::size_t i = 0; i < ndims; ++i) {
      boost::python::slice slice =
        boost::python::extract<boost::python::slice>(selection[i]);
      int slice_start = boost::python::extract<int>(slice.start());
      int slice_stop = boost::python::extract<int>(slice.stop());
      int slice_step = boost::python::extract<int>(slice.step());
      DXTBX_ASSERT(slice_step == 1);
      DXTBX_ASSERT(slice_stop > slice_start);
      dims[i] = slice_stop - slice_start;
      start[i] = slice_start;
      count[i] = dims[i];
      DXTBX_ASSERT(start[i] + count[i] <= dataset_dims[i]);
    }

    // Create the data array
    scitbx::af::flex_grid<> grid(dims);
    scitbx::af::versa<T, scitbx::af::flex_grid<> > data(
      grid, scitbx::af::init_functor_null<T>());

    // Create the dataspace id
    herr_t status1 = H5Sselect_hyperslab(
      file_space_id, H5S_SELECT_SET, &start[0], NULL, &count[0], NULL);
    DXTBX_ASSERT(status1 >= 0);

    // Create the memory space size
    hid_t mem_space_id = H5Screate_simple(ndims, &count[0], NULL);

    // Copy the data
    herr_t status2 = custom_read<T>(dataset_id, mem_space_id, file_space_id, data);
    DXTBX_ASSERT(status2 >= 0);

    // Close some stuff
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);

    // Return the data
    return data;
  }

  BOOST_PYTHON_MODULE(dxtbx_format_nexus_ext) {
    def("dataset_as_flex_int", &dataset_as_flex<int>);
    def("dataset_as_flex_double", &dataset_as_flex<double>);
    def("dataset_as_flex_float", &dataset_as_flex<float>);
    def("image_as_flex_int", &image_as_flex<int>);
    def("image_as_flex_double", &image_as_flex<double>);
    def("image_as_flex_float", &image_as_flex<float>);
  }

}}}  // namespace dxtbx::format::boost_python
