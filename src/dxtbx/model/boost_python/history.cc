#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dxtbx/model/history.h>

namespace dxtbx { namespace model { namespace boost_python {

  void export_history() {
    boost::python::class_<History>("History")
      .def("set_history", &History::set_history_from_list)
      .def("get_history", &History::get_history_as_list);
  }
}}}  // namespace dxtbx::model::boost_python