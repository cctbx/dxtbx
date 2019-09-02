
#include <memory>

#include <stdio.h>

#include <boost/python.hpp>
#include <boost/move/unique_ptr.hpp>

#include <cbf.h>

namespace py = boost::python;

struct PySwigObject {
  PyObject_HEAD void *ptr;
  const char *desc;
};

// Unwrapping SWIG pointer from
// https://wiki.python.org/moin/boost.python/HowTo#SWIG_exposed_C.2B-.2B-_object_from_Python
void *extract_swig_wrapped_pointer(PyObject *obj) {
  char thisStr[] = "this";
  // first we need to get the this attribute from the Python Object
  if (!PyObject_HasAttrString(obj, thisStr)) return NULL;

  PyObject *thisAttr = PyObject_GetAttrString(obj, thisStr);
  if (thisAttr == NULL) return NULL;
  // This Python Object is a SWIG Wrapper and contains our pointer
  void *pointer = ((PySwigObject *)thisAttr)->ptr;
  Py_DECREF(thisAttr);
  return pointer;
}

namespace dxtbx { namespace format { namespace boost_python {

  /// Access the internal buffer and pass it to cbflib
  py::object cbf_read_buffer(py::object handle, py::object data, int flags = 0) {
    if (!PyBytes_Check(data.ptr())) {
      PyErr_SetString(PyExc_ValueError, "buffer must be a bytes-like object");
      py::throw_error_already_set();
    }

    // Extract the opaque CBF object from the SWIG wrapper
    cbf_handle_struct *cbf_handle =
      reinterpret_cast<cbf_handle_struct *>(extract_swig_wrapped_pointer(handle.ptr()));

    int buffer_length = PyBytes_Size(data.ptr());
    char *buffer = PyBytes_AsString(data.ptr());

    int err = cbf_read_buffered_file(
      cbf_handle, NULL /*nullptr*/, flags, buffer, buffer_length);

    if (err) {
      PyErr_Format(PyExc_RuntimeError, "cbflib read_file returned error %d", err);
      py::throw_error_already_set();
    }
    return data;
  }

  /// Declare the cbf-reading classes
  void export_cbf_read_buffer() {
    using namespace boost::python;

    def("cbf_read_buffer",
        cbf_read_buffer,
        args("handle", "data", "flags"),
        "Open a buffer as a CBF file with CBFlib.\n\n"
        "Args:\n"
        "    handle (pycbf.cbf_handle_struct): The CBF handle object\n"
        "    data (bytes): The data buffer containing the CBF file\n"
        "    flags (int): The flags to pass to");
  }
}}}  // namespace dxtbx::format::boost_python
