// Hide a warning in boost on including this
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/python.hpp>
#include <boost/optional.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <cstdint>
#include <complex>
#include <memory>
#include <assert.h>

#include "Python.h"

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/shared_plain.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>

using boost::optional;
using std::cout;
using std::endl;

namespace py = pybind11;
namespace bpy = boost::python;
namespace af = scitbx::af;

template <typename T>
using grid = af::versa<T, af::flex_grid<>>;

// We can do this statically with templates, but makes code somewhat impenetrable
// and TMP-heavy. Just have a macro.
#define ALL_FLEX_TYPES                                                                \
  grid<bool>, grid<size_t>, grid<uint8_t>, grid<uint16_t>, grid<uint32_t>, grid<int>, \
    grid<long>, grid<int8_t>, grid<int16_t>, grid<int32_t>, grid<int64_t>,            \
    grid<float>, grid<double>, grid<scitbx::vec3<double>>, grid<scitbx::vec3<int>>,   \
    grid<scitbx::vec2<double>>, grid<scitbx::mat3<double>>,                           \
    grid<af::tiny<std::size_t, 2>>, grid<std::complex<double>>
// Unwrapped, and possibly doesn't make sense:
// void wrap_flex_std_string(); - differently sized per element
// void wrap_flex_sym_mat3_double(); - nonlinear memory layout

const std::invalid_argument ERR_NON_CONTIGUOUS{
  "numpy array is non-c-contiguous - flex arrays must be c-contiguous"};

bool is_array_c_contiguous(py::array array) {
  return array.attr("flags")["C_CONTIGUOUS"].cast<bool>();
}

template <typename T>
py::buffer_info get_buffer_specific(grid<T> flex) {
  // Build the strides for each dimension by iterating over the sub-dimensions
  std::vector<size_t> strides;
  for (int i = 0; i < flex.accessor().nd(); ++i) {
    auto stride = sizeof(T);
    for (int j = i + 1; j < flex.accessor().nd(); ++j) {
      stride *= flex.accessor().all()[j];
    }
    strides.push_back(stride);
  }
  return py::buffer_info(&flex.front(),
                         sizeof(T),
                         py::format_descriptor<T>::format(),
                         flex.accessor().nd(),
                         flex.accessor().all(),
                         strides);
}

template <typename T>
py::buffer_info get_buffer_specific(grid<scitbx::vec3<T>> flex) {
  std::vector<size_t> dim_sizes;
  for (auto size : flex.accessor().all()) {
    dim_sizes.push_back(size);
  }
  dim_sizes.push_back(3);

  std::vector<size_t> strides;
  for (int i = 0; i < dim_sizes.size(); ++i) {
    auto stride = sizeof(T);
    for (int j = i + 1; j < dim_sizes.size(); ++j) {
      stride *= dim_sizes[j];
    }
    strides.push_back(stride);
  }

  return py::buffer_info(&flex.front(),
                         sizeof(T),
                         py::format_descriptor<T>::format(),
                         dim_sizes.size(),
                         dim_sizes,
                         strides);
}

template <typename T>
py::buffer_info get_buffer_specific(grid<scitbx::vec2<T>> flex) {
  std::vector<size_t> dim_sizes;
  for (auto size : flex.accessor().all()) {
    dim_sizes.push_back(size);
  }
  dim_sizes.push_back(2);

  std::vector<size_t> strides;
  for (int i = 0; i < dim_sizes.size(); ++i) {
    auto stride = sizeof(T);
    for (int j = i + 1; j < dim_sizes.size(); ++j) {
      stride *= dim_sizes[j];
    }
    strides.push_back(stride);
  }

  return py::buffer_info(&flex.front(),
                         sizeof(T),
                         py::format_descriptor<T>::format(),
                         dim_sizes.size(),
                         dim_sizes,
                         strides);
}

template <typename T>
py::buffer_info get_buffer_specific(grid<scitbx::mat3<T>> flex) {
  std::vector<size_t> dim_sizes;
  for (auto size : flex.accessor().all()) {
    dim_sizes.push_back(size);
  }
  // Flex exposes matrix as a 1-d array of nine values - let's expose to
  // numpy as a 3x3 array. We can always get this back with np.reshape(X,9)
  dim_sizes.push_back(3);
  dim_sizes.push_back(3);

  std::vector<size_t> strides;
  for (int i = 0; i < dim_sizes.size(); ++i) {
    auto stride = sizeof(T);
    for (int j = i + 1; j < dim_sizes.size(); ++j) {
      stride *= dim_sizes[j];
    }
    strides.push_back(stride);
  }

  return py::buffer_info(&flex.front(),
                         sizeof(T),
                         py::format_descriptor<T>::format(),
                         dim_sizes.size(),
                         dim_sizes,
                         strides);
}

template <typename T>
py::buffer_info get_buffer_specific(grid<af::tiny<T, 2>> flex) {
  std::vector<size_t> dim_sizes;
  for (auto size : flex.accessor().all()) {
    dim_sizes.push_back(size);
  }
  dim_sizes.push_back(2);

  std::vector<size_t> strides;
  for (int i = 0; i < dim_sizes.size(); ++i) {
    auto stride = sizeof(T);
    for (int j = i + 1; j < dim_sizes.size(); ++j) {
      stride *= dim_sizes[j];
    }
    strides.push_back(stride);
  }

  return py::buffer_info(&flex.front(),
                         sizeof(T),
                         py::format_descriptor<T>::format(),
                         dim_sizes.size(),
                         dim_sizes,
                         strides);
}

template <typename T>
py::buffer_info get_buffer_specific(grid<std::complex<T>> flex) {
  std::vector<size_t> dim_sizes;
  for (auto size : flex.accessor().all()) {
    dim_sizes.push_back(size);
  }

  std::vector<size_t> strides;
  for (int i = 0; i < dim_sizes.size(); ++i) {
    auto stride = sizeof(T) * 2;
    for (int j = i + 1; j < dim_sizes.size(); ++j) {
      stride *= dim_sizes[j];
    }
    strides.push_back(stride);
  }

  return py::buffer_info(&flex.front(),
                         sizeof(T) * 2,
                         "Zd",
                         //  py::format_descriptor<T>::format(),
                         dim_sizes.size(),
                         dim_sizes,
                         strides);
}

/// Get the buffer info for a single flex object type
template <typename T>
py::buffer_info get_buffer(PyObject *ptr) {
  if (bpy::extract<T>(ptr).check()) {
    T obj = bpy::extract<T>(ptr);
    return get_buffer_specific(obj);
  }
  throw std::invalid_argument("Could not determine flex class of object");
}

/// Get the buffer info from multiple flex object types
template <typename T, typename Q, typename... Args>
py::buffer_info get_buffer(PyObject *ptr) {
  if (bpy::extract<T>(ptr).check()) {
    T obj = bpy::extract<T>(ptr);
    return get_buffer_specific(obj);
  }
  return get_buffer<Q, Args...>(ptr);
}

/// Get a buffer_info object for a
py::buffer_info get_buffer(py::object obj) {
  return get_buffer<ALL_FLEX_TYPES>(obj.ptr());
}

template <typename T>
optional<af::sharing_handle *> get_sharing_handle(PyObject *ptr) {
  if (bpy::extract<T>(ptr).check()) {
    T obj = bpy::extract<T>(ptr);
    return obj.handle();
  }
  return boost::none;
}

template <typename T, typename Q, typename... Args>
optional<af::sharing_handle *> get_sharing_handle(PyObject *ptr) {
  if (bpy::extract<T>(ptr).check()) {
    T obj = bpy::extract<T>(ptr);
    return obj.handle();
  }
  return get_sharing_handle<Q, Args...>(ptr);
}

optional<af::sharing_handle *> get_sharing_handle(py::object obj) {
  return get_sharing_handle<ALL_FLEX_TYPES>(obj.ptr());
}
/** Acts as a buffer wrapper for scitbx objects.
 *
 * This holds on to the reference to the object, and provides a buffer
 * interface for passing to numpy, or any other interface that
 * understands PEP-3118-style buffers.
 *
 * */
class Scuffer {
public:
  Scuffer(py::object flex_obj) : _flex_obj(flex_obj) {
    // Call the buffer routine here so that we throw errors early.
    get_buffer_info();
  }

  py::buffer_info get_buffer_info() {
    return get_buffer(_flex_obj);
  }

  py::object base() {
    return _flex_obj;
  }

private:
  py::object _flex_obj;
};

/** A array_family sharing_handle that owns the memory via a numpy base
 * **/
class numpy_sharing_handle : public af::sharing_handle {
public:
  numpy_sharing_handle(py::array np_array) : _np_array(np_array) {
    size = capacity = _np_array.nbytes();
    data = reinterpret_cast<char *>(_np_array.mutable_data());
  }

  virtual ~numpy_sharing_handle() {
    deallocate();
  }
  virtual void deallocate() {
    // Since we don't own the memory,
    _np_array = py::none();
    size = capacity = 0;
    data = nullptr;
  }

  virtual void swap(sharing_handle &other) {
    auto other_np = dynamic_cast<numpy_sharing_handle *>(&other);
    if (other_np == nullptr)
      throw std::invalid_argument("Cannot swap numpy and non-numpy handles");
    std::swap(other_np->_np_array, _np_array);
    af::sharing_handle::swap(other);
  }

  virtual py::array base() {
    return _np_array;
  }

private:
  py::array _np_array;
};

/** Convert a flex object directly to numpy
 *
 * Without creating an intermediate object. This creates a Scuffer
 * object under-the-hood and passes it to numpy.
 *
 * */
py::object to_numpy(py::object flex_array) {
  // If we've been given a numpy array already... just do nothing.
  if (py::isinstance<py::array>(flex_array)) {
    return flex_array;
  }
  // We want to try and work out if this flex object is itself wrapping
  // numpy, because then we can just return the underlying numpy object.
  if (auto sh = get_sharing_handle(flex_array)) {
    if (auto as_numpy_handle = dynamic_cast<numpy_sharing_handle *>(*sh)) {
      return as_numpy_handle->base();
    }
  }
  py::object scuffer = py::cast(Scuffer(flex_array));
  auto numpy = py::module_::import("numpy");
  return numpy.attr("asarray")(scuffer);
}

/** Convert a numpy object to a typed array_family object.
 *
 * Args:
 *  np_array: The numpy python object
 *  ignore_dims: The number of dimensions to ignore. If set, this means
 *    that the target type itself contains the objects (e.g. vec3 or mat3)
 *    and so they don't need to be sized into the flex_grid.
 * **/
template <typename T>
py::object numpy_to_array_family(py::array np_array, int ignore_dims = 0) {
  auto handle = new numpy_sharing_handle(np_array);
  if (np_array.ndim() > 10 + ignore_dims) {
    throw std::invalid_argument("Default flex grid only supports up to 10 dimensions");
  }
  assert(ignore_dims < np_array.ndim());

  // Work out the grid size
  auto dims = typename T::accessor_type::index_type();
  for (int i = 0; i < np_array.ndim() - ignore_dims; ++i) {
    dims.push_back(np_array.shape(i));
  }
  auto grid = typename T::accessor_type(dims);
  auto scitbx_obj = T(handle, grid);
  auto boost_obj = bpy::object(scitbx_obj);
  // Release our interest in the shared handle
  handle->use_count--;
  return py::reinterpret_borrow<py::object>(boost_obj.ptr());
}

py::object from_numpy(py::object array) {
  // Firstly, check if we are wrapping a flex object already
  // now, let's check that this isn't already a flex object
  try {
    get_buffer(array);
    // This isn't a numpy array - it's already a known flex object!
    return array;
  } catch (const std::invalid_argument &) {
    // This isn't a flex object(directly)
  }
  if (!py::isinstance<py::array>(array)) {
    throw std::invalid_argument(
      "Cannot currently convert from non-numpy array format to flex");
  }
  auto np_array = py::array(array);

  // Check that this array is contiguous
  if (!is_array_c_contiguous(np_array)) {
    throw ERR_NON_CONTIGUOUS;
  }

  // If this was directly converted from a flex array, give back the original
  // object; In any other case we want to wrap, because of slicing/metadata
  if (np_array.base()) {
    if (py::isinstance<py::memoryview>(np_array.base())) {
      if (py::isinstance<Scuffer>(np_array.base().attr("obj"))) {
        // Ah, this came from flex originally
        auto scuffer = np_array.base().attr("obj").cast<Scuffer &>();
        return scuffer.base();
      }
    }
  }

  // Check that we recognise this type
  std::string known_types = "BHILQbhilqdDf?";
  auto dtype_obj = np_array.attr("dtype");
  auto dtype = dtype_obj.attr("char").cast<char>();
  if (known_types.find(dtype) == std::string::npos) {
    throw std::invalid_argument(std::string("Unrecognised numpy array dtype '") + dtype
                                + "'");
  }

  // Need to resolve mappings - some of flex is bound as specific integer
  // sizing, some is bound directly to the sized type. We need to resolve
  // the differences when we've been given a numpy dtype that doesn't
  // have an exact corresponding flex implementation.
  //
  // Flex binds only integer types:
  //    bool, uint8, uint16, uint32, size_t, int, long, int8, int16, [int64 - win]
  // so the direct conversions are:
  //    ? -> flex.bool
  //    i -> flex.int
  //    l -> flex.long
  //    f -> flex.float
  //    d -> flex.double
  //    D -> flex.complex
  // everything else needs to be mapped via size/signededness; although size_t
  // is mapped directly, it doesn't have a fixed numpy format character.

  // Firstly, try to map the direct equivalents
  if (dtype == '?') {
    return numpy_to_array_family<af::versa<bool, af::flex_grid<>>>(np_array);
  } else if (dtype == 'i') {
    return numpy_to_array_family<af::versa<int, af::flex_grid<>>>(np_array);
  } else if (dtype == 'l' && (sizeof(long) == sizeof(int))) {
    // In cases where int and long are degenerate (windows), prefer int
    return numpy_to_array_family<af::versa<int, af::flex_grid<>>>(np_array);
  } else if (dtype == 'l') {
    return numpy_to_array_family<af::versa<long, af::flex_grid<>>>(np_array);
  } else if (dtype == 'f') {
    return numpy_to_array_family<af::versa<float, af::flex_grid<>>>(np_array);
  } else if (dtype == 'd') {
    return numpy_to_array_family<af::versa<double, af::flex_grid<>>>(np_array);
  } else if (dtype == 'D') {
    return numpy_to_array_family<af::versa<std::complex<double>, af::flex_grid<>>>(
      np_array);
  }

  // Extract information about the dtype so we can match it up
  auto numpy = py::module_::import("numpy");
  auto is_integer =
    numpy.attr("issubdtype")(dtype_obj, numpy.attr("integer")).cast<bool>();
  // We _think_ that we've covered this, but double check
  if (!is_integer) {
    throw std::runtime_error(std::string("Unknown/unhandled non-integer data type: ")
                             + dtype);
  }

  auto is_signed =
    numpy.attr("issubdtype")(dtype_obj, numpy.attr("signedinteger")).cast<bool>();
  auto itemsize = dtype_obj.attr("itemsize").cast<int>();

  // Now, work out which data type to use based on metadata about the type
  if (is_signed) {
    if (itemsize == sizeof(int)) {
      return numpy_to_array_family<af::versa<int, af::flex_grid<>>>(np_array);
    } else if (itemsize == sizeof(long)) {
      return numpy_to_array_family<af::versa<long, af::flex_grid<>>>(np_array);
    } else if (itemsize == sizeof(int8_t)) {
      return numpy_to_array_family<af::versa<int8_t, af::flex_grid<>>>(np_array);
    } else if (itemsize == sizeof(int16_t)) {
      return numpy_to_array_family<af::versa<int16_t, af::flex_grid<>>>(np_array);
    } else if (itemsize == sizeof(int64_t)) {
      // Only hit this in case where long != uint64, and so this should be bound
      return numpy_to_array_family<af::versa<int64_t, af::flex_grid<>>>(np_array);
    }
  } else {
    // Unsigned -
    if (itemsize == sizeof(size_t)) {
      return numpy_to_array_family<af::versa<size_t, af::flex_grid<>>>(np_array);
    } else if (itemsize == sizeof(uint8_t)) {
      return numpy_to_array_family<af::versa<uint8_t, af::flex_grid<>>>(np_array);
    } else if (itemsize == sizeof(uint16_t)) {
      return numpy_to_array_family<af::versa<uint16_t, af::flex_grid<>>>(np_array);
    } else if (itemsize == sizeof(uint32_t)) {
      return numpy_to_array_family<af::versa<uint32_t, af::flex_grid<>>>(np_array);
    }
  }
  throw std::runtime_error(std::string("Unhandled integer array type '") + dtype
                           + "' - unrecognised size " + std::to_string(itemsize));
}

/// More structured arrays need to be explicitly requested
template <template <class> class VecType>
py::object vec_from_numpy(py::array np_array) {
  static_assert(VecType<int>::fixed_size == 2 || VecType<int>::fixed_size == 3,
                "Only vec2/vec3 supported");
  // Only accept arrays whose last dimension is the size of this object
  if (np_array.shape(np_array.ndim() - 1) != VecType<int>::fixed_size) {
    throw std::invalid_argument("Input array last dimension is not size "
                                + std::to_string(VecType<int>::fixed_size));
  }

  auto dtype = np_array.attr("dtype").attr("char").cast<char>();

  std::string accepted_types = VecType<int>::fixed_size == 2 ? "dQ" : "di";
  if (accepted_types.find(dtype) == std::string::npos) {
    throw std::invalid_argument(std::string("Can only convert ") + accepted_types
                                + std::string(" to vec - not '") + dtype + "'");
  }
  if (dtype == 'd') {
    return numpy_to_array_family<af::versa<VecType<double>, af::flex_grid<>>>(np_array,
                                                                              1);
  } else if (dtype == 'i') {
    return numpy_to_array_family<af::versa<VecType<int>, af::flex_grid<>>>(np_array, 1);
  } else if (dtype == 'Q') {
    // The only known case here is flex.tiny_size_t_2, which is bound as
    //    versa<tiny<size_t,2> >
    // for some reason instead of
    //    versa<vec2<size_t> >
    // See the Zen: Consistency. Rather than working out how to contort
    // the template system into bending over backwards for this case that
    // we probably won't ever use, override the template parameter and
    // just directly bind it. I'm very sorry if this has cause a headache
    // along the lines of "Why is it giving me a tiny<> when I ask for a
    // vec<>!?!?!?".
    return numpy_to_array_family<af::versa<af::tiny<size_t, 2>, af::flex_grid<>>>(
      np_array, 1);
  }
  throw std::runtime_error(std::string("Unhandled array type '") + dtype + "'");
}

/// Decide which sized vector we want to convert to, and hand off to the
/// specialization
py::object vecs_from_numpy(py::array np_array) {
  // Check that this array is contiguous
  if (!is_array_c_contiguous(np_array)) {
    throw ERR_NON_CONTIGUOUS;
  }

  if (np_array.shape(np_array.ndim() - 1) == 3) {
    return vec_from_numpy<scitbx::vec3>(np_array);
  } else if (np_array.shape(np_array.ndim() - 1) == 2) {
    return vec_from_numpy<scitbx::vec2>(np_array);
  }
  throw std::invalid_argument(
    "Invalid input array: last numpy dimension must be 2 or 3 to convert to vector");
}

py::object mat3_from_numpy(py::array np_array) {
  // Check that this array is contiguous
  if (!is_array_c_contiguous(np_array)) {
    throw ERR_NON_CONTIGUOUS;
  }

  auto nd = np_array.ndim();
  // Check our last dimension(s) are either x9 or x3x3
  bool last_is_9 = np_array.shape(nd - 1) == 9;
  bool last_is_3x3 =
    nd > 2 && np_array.shape(nd - 1) == 3 && np_array.shape(nd - 2) == 3;
  if (!(last_is_3x3 || last_is_9)) {
    throw std::invalid_argument(
      "Input array is not ...x3x3 or ...x9 - can't derive target shape");
  }
  auto dtype = np_array.attr("dtype").attr("char").cast<char>();
  if (dtype != 'd') {
    throw std::invalid_argument(
      std::string("Only mat3_double is known bound in flex - cannot convert '")
      + std::to_string(dtype) + "'");
  }
  if (last_is_9) {
    return numpy_to_array_family<af::versa<scitbx::mat3<double>, af::flex_grid<>>>(
      np_array, 1);
  } else if (last_is_3x3) {
    return numpy_to_array_family<af::versa<scitbx::mat3<double>, af::flex_grid<>>>(
      np_array, 2);
  }
  throw std::runtime_error("Unhandled matrix size");
}

PYBIND11_MODULE(dxtbx_flumpy, m) {
  py::class_<Scuffer>(m, "Scuffer", py::buffer_protocol())
    .def(py::init<py::object>())
    .def_property("base", &Scuffer::base, nullptr)
    .def_buffer(&Scuffer::get_buffer_info);
  m.def("to_numpy",
        &to_numpy,
        "Convert a flex object into a numpy array with zero copying");
  m.def("from_numpy",
        &from_numpy,
        "Convert a numpy object into the equivalent (flat) flex array");
  m.def("vec_from_numpy",
        &vecs_from_numpy,
        "Convert a numpy object to a flex.vec2 or .vec3, depending on input array");
  m.def("mat3_from_numpy", &mat3_from_numpy, "Convert a numpy object to a flex.mat3");

  // Make sure that we have imported flex - cannot do boost::python conversions
  // otherwise
  pybind11::module::import("scitbx.array_family.flex");
}
