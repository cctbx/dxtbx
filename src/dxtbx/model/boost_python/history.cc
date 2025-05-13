#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dxtbx/model/history.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model { namespace boost_python {

  struct HistoryPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getstate(boost::python::object obj) {
      const History &history = boost::python::extract<const History &>(obj)();
      return boost::python::make_tuple(obj.attr("__dict__"),
                                       history.get_history_as_list());
    }

    static void setstate(boost::python::object obj, boost::python::tuple state) {
      History &history = boost::python::extract<History &>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 2);

      // restore the object's __dict__
      boost::python::dict d =
        boost::python::extract<boost::python::dict>(obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      history.set_history_from_list(
        boost::python::extract<boost::python::list>(state[1])());
    }
  };

  void export_history() {
    boost::python::class_<History>("History")
      .def("set_history", &History::set_history_from_list)
      .def("get_history", &History::get_history_as_list)
      .def(
        "append_history_item",
        &History::append_history_item(arg("dispatcher"), arg("version"), arg("flag")))
      .def_pickle(HistoryPickleSuite());
  }
}}}  // namespace dxtbx::model::boost_python