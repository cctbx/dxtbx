/*
 * experiment_list.cc
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/slice.hpp>
#include <string>
#include <sstream>
#include <scitbx/constants.h>
#include <scitbx/boost_python/slice.h>
#include <scitbx/boost_python/utils.h>
#include <dxtbx/model/experiment_list.h>
#include <dxtbx/model/boost_python/to_from_dict.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  template<typename Beam>
  static ExperimentList<Beam> *make_experiment_list(boost::python::object items) {
    ExperimentList<Beam> *experiment = new ExperimentList<Beam>();

    for (std::size_t i = 0; i < boost::python::len(items); ++i) {
      experiment->append(boost::python::extract<Experiment<Beam>>(items[i])());
    }

    return experiment;
  }

  template<typename Beam>
  struct ExperimentListPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(const ExperimentList<Beam> &obj) {
      boost::python::list experiments;
      for (std::size_t i = 0; i < obj.size(); ++i) {
        experiments.append(obj[i]);
      }
      return boost::python::make_tuple(experiments);
    }
  };

  /**
   * Return function pointers to overrides for different types
   */
  template<typename Beam>
  struct experiment_list_contains_pointers {
    typedef bool (ExperimentList<Beam>::*beam_type)(
      const boost::shared_ptr<Beam> &) const;
    typedef bool (ExperimentList<Beam>::*detector_type)(
      const boost::shared_ptr<Detector> &) const;
    typedef bool (ExperimentList<Beam>::*goniometer_type)(
      const boost::shared_ptr<Goniometer> &) const;
    typedef bool (ExperimentList<Beam>::*scan_type)(const boost::shared_ptr<Scan> &) const;
    typedef bool (ExperimentList<Beam>::*crystal_type)(
      const boost::shared_ptr<CrystalBase> &) const;
    typedef bool (ExperimentList<Beam>::*object_type)(boost::python::object) const;

    static beam_type beam() {
      return &ExperimentList<Beam>::contains;
    }

    static detector_type detector() {
      return &ExperimentList<Beam>::contains;
    }

    static goniometer_type goniometer() {
      return &ExperimentList<Beam>::contains;
    }

    static scan_type scan() {
      return &ExperimentList<Beam>::contains;
    }

    static crystal_type crystal() {
      return &ExperimentList<Beam>::contains;
    }

    static object_type object() {
      return &ExperimentList<Beam>::contains;
    }
  };

  /**
   * Return function pointers to overrides for different types
   */
  template<typename Beam>
  struct experiment_list_replace_pointers {
    typedef void (ExperimentList<Beam>::*beam_type)(boost::shared_ptr<Beam>,
                                              boost::shared_ptr<Beam>);
    typedef void (ExperimentList<Beam>::*detector_type)(boost::shared_ptr<Detector>,
                                                  boost::shared_ptr<Detector>);
    typedef void (ExperimentList<Beam>::*goniometer_type)(boost::shared_ptr<Goniometer>,
                                                    boost::shared_ptr<Goniometer>);
    typedef void (ExperimentList<Beam>::*scan_type)(boost::shared_ptr<Scan>,
                                              boost::shared_ptr<Scan>);
    typedef void (ExperimentList<Beam>::*crystal_type)(boost::shared_ptr<CrystalBase>,
                                                 boost::shared_ptr<CrystalBase>);
    typedef void (ExperimentList<Beam>::*object_type)(boost::python::object,
                                                boost::python::object);

    static beam_type beam() {
      return &ExperimentList<Beam>::replace;
    }

    static detector_type detector() {
      return &ExperimentList<Beam>::replace;
    }

    static goniometer_type goniometer() {
      return &ExperimentList<Beam>::replace;
    }

    static scan_type scan() {
      return &ExperimentList<Beam>::replace;
    }

    static crystal_type crystal() {
      return &ExperimentList<Beam>::replace;
    }

    static object_type object() {
      return &ExperimentList<Beam>::replace;
    }
  };

  /**
   * Return function pointers to overrides for different types
   */
  template<typename Beam>
  struct experiment_list_indices_pointers {
    typedef scitbx::af::shared<std::size_t> (ExperimentList<Beam>::*beam_type)(
      const boost::shared_ptr<Beam> &) const;
    typedef scitbx::af::shared<std::size_t> (ExperimentList<Beam>::*detector_type)(
      const boost::shared_ptr<Detector> &) const;
    typedef scitbx::af::shared<std::size_t> (ExperimentList<Beam>::*goniometer_type)(
      const boost::shared_ptr<Goniometer> &) const;
    typedef scitbx::af::shared<std::size_t> (ExperimentList<Beam>::*scan_type)(
      const boost::shared_ptr<Scan> &) const;
    typedef scitbx::af::shared<std::size_t> (ExperimentList<Beam>::*crystal_type)(
      const boost::shared_ptr<CrystalBase> &) const;
    typedef scitbx::af::shared<std::size_t> (ExperimentList<Beam>::*object_type)(
      boost::python::object) const;

    static beam_type beam() {
      return &ExperimentList<Beam>::indices;
    }

    static detector_type detector() {
      return &ExperimentList<Beam>::indices;
    }

    static goniometer_type goniometer() {
      return &ExperimentList<Beam>::indices;
    }

    static scan_type scan() {
      return &ExperimentList<Beam>::indices;
    }

    static crystal_type crystal() {
      return &ExperimentList<Beam>::indices;
    }

    static object_type object() {
      return &ExperimentList<Beam>::indices;
    }
  };

  template<typename Beam>
  void experiment_list_setitem(ExperimentList<Beam> &self, int n, Experiment<Beam> item) {
    if (n < 0) {
      n += self.size();
    }
    if (n >= self.size() || n < 0) {
      scitbx::boost_python::raise_index_error();
    }
    self[n] = item;
  }

  template<typename Beam>
  Experiment<Beam> &experiment_list_getitem(ExperimentList<Beam> &self, int n) {
    if (n < 0) {
      n += self.size();
    }
    if (n >= self.size() || n < 0) {
      scitbx::boost_python::raise_index_error();
    }
    return self[n];
  }

  template<typename Beam>
  ExperimentList<Beam> experiment_list_getitem_slice(const ExperimentList<Beam> &self, slice s) {
    scitbx::boost_python::adapted_slice as(s, self.size());
    ExperimentList<Beam> result;
    for (std::size_t i = as.start; i < as.stop && i < self.size(); i += as.step) {
      result.append(self[i]);
    }
    return result;
  }

  template<typename Beam>
  void experiment_list_delitem(ExperimentList<Beam> &self, int n) {
    if (n < 0) {
      n += self.size();
    }
    if (n >= self.size() || n < 0) {
      scitbx::boost_python::raise_index_error();
    }
    self.erase(n);
  }

  void export_experiment_list() {
    class_<ExperimentList<MonochromaticBeam> >("ExperimentList")
      .def("__init__", make_constructor(&make_experiment_list<MonochromaticBeam>, default_call_policies()))
      .def("identifiers", &ExperimentList<MonochromaticBeam>::identifiers)
      .def("find", &ExperimentList<MonochromaticBeam>::find)
      .def("append", &ExperimentList<MonochromaticBeam>::append)
      .def("extend", &ExperimentList<MonochromaticBeam>::extend)
      .def("clear", &ExperimentList<MonochromaticBeam>::clear)
      .def("empty", &ExperimentList<MonochromaticBeam>::empty)
      .def("__getitem__", &experiment_list_getitem<MonochromaticBeam>, return_internal_reference<>())
      .def("__getitem__", &experiment_list_getitem_slice<MonochromaticBeam>)
      .def("__setitem__", &experiment_list_setitem<MonochromaticBeam>)
      .def("__delitem__", &experiment_list_delitem<MonochromaticBeam>)
      .def("__iter__", iterator<ExperimentList<MonochromaticBeam>, return_internal_reference<> >())
      .def("__contains__", experiment_list_contains_pointers<MonochromaticBeam>::beam())
      .def("__contains__", experiment_list_contains_pointers<MonochromaticBeam>::detector())
      .def("__contains__", experiment_list_contains_pointers<MonochromaticBeam>::goniometer())
      .def("__contains__", experiment_list_contains_pointers<MonochromaticBeam>::scan())
      .def("__contains__", experiment_list_contains_pointers<MonochromaticBeam>::crystal())
      .def("__contains__", experiment_list_contains_pointers<MonochromaticBeam>::object())
      .def("replace", experiment_list_replace_pointers<MonochromaticBeam>::beam())
      .def("replace", experiment_list_replace_pointers<MonochromaticBeam>::detector())
      .def("replace", experiment_list_replace_pointers<MonochromaticBeam>::goniometer())
      .def("replace", experiment_list_replace_pointers<MonochromaticBeam>::scan())
      .def("replace", experiment_list_replace_pointers<MonochromaticBeam>::crystal())
      .def("replace", experiment_list_replace_pointers<MonochromaticBeam>::object())
      .def("indices", experiment_list_indices_pointers<MonochromaticBeam>::beam())
      .def("indices", experiment_list_indices_pointers<MonochromaticBeam>::detector())
      .def("indices", experiment_list_indices_pointers<MonochromaticBeam>::goniometer())
      .def("indices", experiment_list_indices_pointers<MonochromaticBeam>::scan())
      .def("indices", experiment_list_indices_pointers<MonochromaticBeam>::crystal())
      .def("indices", experiment_list_indices_pointers<MonochromaticBeam>::object())
      .def("remove_on_experiment_identifiers",
           &ExperimentList<MonochromaticBeam>::remove_on_experiment_identifiers)
      .def("select_on_experiment_identifiers",
           &ExperimentList<MonochromaticBeam>::select_on_experiment_identifiers)
      .def("where",
           &ExperimentList<MonochromaticBeam>::where,
           (arg("beam") = boost::shared_ptr<MonochromaticBeam>(),
            arg("detector") = boost::shared_ptr<Detector>(),
            arg("goniometer") = boost::shared_ptr<Goniometer>(),
            arg("scan") = boost::shared_ptr<Scan>(),
            arg("crystal") = boost::shared_ptr<CrystalBase>(),
            arg("profile") = boost::python::object(),
            arg("imageset") = boost::python::object(),
            arg("scaling_model") = boost::python::object()))
      .def("is_consistent", &ExperimentList<MonochromaticBeam>::is_consistent)
      .def("__len__", &ExperimentList<MonochromaticBeam>::size)
      .def_pickle(ExperimentListPickleSuite<MonochromaticBeam>());
    class_<ExperimentList<TOFBeam> >("ExperimentList")
      .def("__init__", make_constructor(&make_experiment_list<TOFBeam>, default_call_policies()))
      .def("identifiers", &ExperimentList<TOFBeam>::identifiers)
      .def("find", &ExperimentList<TOFBeam>::find)
      .def("append", &ExperimentList<TOFBeam>::append)
      .def("extend", &ExperimentList<TOFBeam>::extend)
      .def("clear", &ExperimentList<TOFBeam>::clear)
      .def("empty", &ExperimentList<TOFBeam>::empty)
      .def("__getitem__", &experiment_list_getitem<TOFBeam>, return_internal_reference<>())
      .def("__getitem__", &experiment_list_getitem_slice<TOFBeam>)
      .def("__setitem__", &experiment_list_setitem<TOFBeam>)
      .def("__delitem__", &experiment_list_delitem<TOFBeam>)
      .def("__iter__", iterator<ExperimentList<TOFBeam>, return_internal_reference<> >())
      .def("__contains__", experiment_list_contains_pointers<TOFBeam>::beam())
      .def("__contains__", experiment_list_contains_pointers<TOFBeam>::detector())
      .def("__contains__", experiment_list_contains_pointers<TOFBeam>::goniometer())
      .def("__contains__", experiment_list_contains_pointers<TOFBeam>::scan())
      .def("__contains__", experiment_list_contains_pointers<TOFBeam>::crystal())
      .def("__contains__", experiment_list_contains_pointers<TOFBeam>::object())
      .def("replace", experiment_list_replace_pointers<TOFBeam>::beam())
      .def("replace", experiment_list_replace_pointers<TOFBeam>::detector())
      .def("replace", experiment_list_replace_pointers<TOFBeam>::goniometer())
      .def("replace", experiment_list_replace_pointers<TOFBeam>::scan())
      .def("replace", experiment_list_replace_pointers<TOFBeam>::crystal())
      .def("replace", experiment_list_replace_pointers<TOFBeam>::object())
      .def("indices", experiment_list_indices_pointers<TOFBeam>::beam())
      .def("indices", experiment_list_indices_pointers<TOFBeam>::detector())
      .def("indices", experiment_list_indices_pointers<TOFBeam>::goniometer())
      .def("indices", experiment_list_indices_pointers<TOFBeam>::scan())
      .def("indices", experiment_list_indices_pointers<TOFBeam>::crystal())
      .def("indices", experiment_list_indices_pointers<TOFBeam>::object())
      .def("remove_on_experiment_identifiers",
           &ExperimentList<TOFBeam>::remove_on_experiment_identifiers)
      .def("select_on_experiment_identifiers",
           &ExperimentList<TOFBeam>::select_on_experiment_identifiers)
      .def("where",
           &ExperimentList<TOFBeam>::where,
           (arg("beam") = boost::shared_ptr<TOFBeam>(),
            arg("detector") = boost::shared_ptr<Detector>(),
            arg("goniometer") = boost::shared_ptr<Goniometer>(),
            arg("scan") = boost::shared_ptr<Scan>(),
            arg("crystal") = boost::shared_ptr<CrystalBase>(),
            arg("profile") = boost::python::object(),
            arg("imageset") = boost::python::object(),
            arg("scaling_model") = boost::python::object()))
      .def("is_consistent", &ExperimentList<TOFBeam>::is_consistent)
      .def("__len__", &ExperimentList<TOFBeam>::size)
      .def_pickle(ExperimentListPickleSuite<TOFBeam>());
  }

}}}  // namespace dxtbx::model::boost_python
