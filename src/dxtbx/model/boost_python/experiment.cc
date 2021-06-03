/*
 * experiment.cc
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
#include <string>
#include <sstream>
#include <scitbx/constants.h>
#include <dxtbx/model/experiment.h>
#include <dxtbx/model/boost_python/to_from_dict.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  template<typename Beam>
  struct ExperimentPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(const Experiment<Beam> &obj) {
      return boost::python::make_tuple(obj.get_beam(),
                                       obj.get_detector(),
                                       obj.get_goniometer(),
                                       obj.get_scan(),
                                       obj.get_crystal(),
                                       obj.get_profile(),
                                       obj.get_imageset(),
                                       obj.get_scaling_model(),
                                       obj.get_identifier());
    }
  };

  /**
   * Return function pointers to overrides for different types
   */
  template<typename Beam>
  struct experiment_contains_pointers {
    typedef bool (Experiment<Beam>::*beam_type)(const boost::shared_ptr<Beam> &) const;
    typedef bool (Experiment<Beam>::*detector_type)(
      const boost::shared_ptr<Detector> &) const;
    typedef bool (Experiment<Beam>::*goniometer_type)(
      const boost::shared_ptr<Goniometer> &) const;
    typedef bool (Experiment<Beam>::*scan_type)(const boost::shared_ptr<Scan> &) const;
    typedef bool (Experiment<Beam>::*crystal_type)(
      const boost::shared_ptr<CrystalBase> &) const;
    typedef bool (Experiment<Beam>::*object_type)(boost::python::object) const;

    static beam_type beam() {
      return &Experiment<Beam>::contains;
    }

    static detector_type detector() {
      return &Experiment<Beam>::contains;
    }

    static goniometer_type goniometer() {
      return &Experiment<Beam>::contains;
    }

    static scan_type scan() {
      return &Experiment<Beam>::contains;
    }

    static crystal_type crystal() {
      return &Experiment<Beam>::contains;
    }

    static object_type object() {
      return &Experiment<Beam>::contains;
    }
  };

  template<typename Beam>
  void export_experiment() {
    class_<Experiment<Beam> >("Experiment")
      .def(init<boost::shared_ptr<Beam>,
                boost::shared_ptr<Detector>,
                boost::shared_ptr<Goniometer>,
                boost::shared_ptr<Scan>,
                boost::shared_ptr<CrystalBase>,
                boost::python::object,
                boost::python::object,
                boost::python::object,
                std::string>((arg("beam") = boost::shared_ptr<Beam>(),
                              arg("detector") = boost::shared_ptr<Detector>(),
                              arg("goniometer") = boost::shared_ptr<Goniometer>(),
                              arg("scan") = boost::shared_ptr<Scan>(),
                              arg("crystal") = boost::shared_ptr<CrystalBase>(),
                              arg("profile") = boost::python::object(),
                              arg("imageset") = boost::python::object(),
                              arg("scaling_model") = boost::python::object(),
                              arg("identifier") = "")))
      .add_property("beam", &Experiment<Beam>::get_beam, &Experiment<Beam>::set_beam)
      .add_property("detector", &Experiment<Beam>::get_detector, &Experiment<Beam>::set_detector)
      .add_property(
        "goniometer", &Experiment<Beam>::get_goniometer, &Experiment<Beam>::set_goniometer)
      .add_property("scan", &Experiment<Beam>::get_scan, &Experiment<Beam>::set_scan)
      .add_property("crystal", &Experiment<Beam>::get_crystal, &Experiment<Beam>::set_crystal)
      .add_property("profile", &Experiment<Beam>::get_profile, &Experiment<Beam>::set_profile)
      .add_property("imageset", &Experiment<Beam>::get_imageset, &Experiment<Beam>::set_imageset)
      .add_property(
        "scaling_model", &Experiment<Beam>::get_scaling_model, &Experiment<Beam>::set_scaling_model)
      .add_property(
        "identifier", &Experiment<Beam>::get_identifier, &Experiment<Beam>::set_identifier)
      .def("__contains__", experiment_contains_pointers<Beam>::beam())
      .def("__contains__", experiment_contains_pointers<Beam>::detector())
      .def("__contains__", experiment_contains_pointers<Beam>::goniometer())
      .def("__contains__", experiment_contains_pointers<Beam>::scan())
      .def("__contains__", experiment_contains_pointers<Beam>::crystal())
      .def("__contains__", experiment_contains_pointers<Beam>::object())
      .def("__eq__", &Experiment<Beam>::operator==)
      .def("is_consistent", &Experiment<Beam>::is_consistent)
      .def("is_still",
           &Experiment<Beam>::is_still,
           "Check if this experiment represents a still image")
      .def("is_sequence",
           &Experiment<Beam>::is_sequence,
           "Check if this experiment represents swept rotation image(s)")
      .def_pickle(ExperimentPickleSuite<Beam>());
  }

}}}  // namespace dxtbx::model::boost_python
