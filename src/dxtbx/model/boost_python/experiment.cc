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
#include <dxtbx/model/beam.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  
  struct ExperimentPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(const Experiment &obj) {
      return boost::python::make_tuple(obj.get_beam(),
                                       obj.get_detector(),
                                       obj.get_goniometer(),
                                       obj.get_sequence(),
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
  
  struct experiment_contains_pointers {
    typedef bool (Experiment::*beam_type)(const boost::python::object &) const;
    typedef bool (Experiment::*detector_type)(
      const boost::shared_ptr<Detector> &) const;
    typedef bool (Experiment::*goniometer_type)(
      const boost::shared_ptr<Goniometer> &) const;
    typedef bool (Experiment::*sequence_type)(const boost::python::object &) const;
    typedef bool (Experiment::*crystal_type)(
      const boost::shared_ptr<CrystalBase> &) const;
    typedef bool (Experiment::*object_type)(boost::python::object) const;

    static beam_type beam() {
      return &Experiment::contains_beam;
    }

    static detector_type detector() {
      return &Experiment::contains;
    }

    static goniometer_type goniometer() {
      return &Experiment::contains;
    }

    static sequence_type sequence() {
      return &Experiment::contains_sequence;
    }

    static crystal_type crystal() {
      return &Experiment::contains;
    }

    static object_type object() {
      return &Experiment::contains;
    }
  };

  void export_experiment() {
    class_<Experiment>("Experiment")
      .def(init<boost::python::object,
                boost::shared_ptr<Detector>,
                boost::shared_ptr<Goniometer>,
                boost::python::object,
                boost::shared_ptr<CrystalBase>,
                boost::python::object,
                boost::python::object,
                boost::python::object,
                std::string>((arg("beam") = boost::python::object(),
                              arg("detector") = boost::shared_ptr<Detector>(),
                              arg("goniometer") = boost::shared_ptr<Goniometer>(),
                              arg("sequence") = boost::python::object(),
                              arg("crystal") = boost::shared_ptr<CrystalBase>(),
                              arg("profile") = boost::python::object(),
                              arg("imageset") = boost::python::object(),
                              arg("scaling_model") = boost::python::object(),
                              arg("identifier") = "")))
      .add_property("beam", &Experiment::get_beam, &Experiment::set_beam)
      .add_property("detector", &Experiment::get_detector, &Experiment::set_detector)
      .add_property(
        "goniometer", &Experiment::get_goniometer, &Experiment::set_goniometer)
      .add_property("sequence", &Experiment::get_sequence, &Experiment::set_sequence)
      .add_property("crystal", &Experiment::get_crystal, &Experiment::set_crystal)
      .add_property("profile", &Experiment::get_profile, &Experiment::set_profile)
      .add_property("imageset", &Experiment::get_imageset, &Experiment::set_imageset)
      .add_property(
        "scaling_model", &Experiment::get_scaling_model, &Experiment::set_scaling_model)
      .add_property(
        "identifier", &Experiment::get_identifier, &Experiment::set_identifier)
      .def("__contains__", experiment_contains_pointers::beam())
      .def("__contains__", experiment_contains_pointers::detector())
      .def("__contains__", experiment_contains_pointers::goniometer())
      .def("__contains__", experiment_contains_pointers::sequence())
      .def("__contains__", experiment_contains_pointers::crystal())
      .def("__contains__", experiment_contains_pointers::object())
      .def("__eq__", &Experiment::operator==)
      .def("is_consistent", &Experiment::is_consistent)
      .def("is_still",
           &Experiment::is_still,
           "Check if this experiment represents a still image")
      .def("is_sequence",
           &Experiment::is_sequence,
           "Check if this experiment represents swept rotation image(s)")
      .def_pickle(ExperimentPickleSuite());
  }

}}}  // namespace dxtbx::model::boost_python
