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

  void export_experiment() {
    class_<Experiment<MonochromaticBeam> >("Experiment")
      .def(init<boost::shared_ptr<MonochromaticBeam>,
                boost::shared_ptr<Detector>,
                boost::shared_ptr<Goniometer>,
                boost::shared_ptr<Scan>,
                boost::shared_ptr<CrystalBase>,
                boost::python::object,
                boost::python::object,
                boost::python::object,
                std::string>((arg("beam") = boost::shared_ptr<MonochromaticBeam>(),
                              arg("detector") = boost::shared_ptr<Detector>(),
                              arg("goniometer") = boost::shared_ptr<Goniometer>(),
                              arg("scan") = boost::shared_ptr<Scan>(),
                              arg("crystal") = boost::shared_ptr<CrystalBase>(),
                              arg("profile") = boost::python::object(),
                              arg("imageset") = boost::python::object(),
                              arg("scaling_model") = boost::python::object(),
                              arg("identifier") = "")))
      .add_property("beam", &Experiment<MonochromaticBeam>::get_beam, &Experiment<MonochromaticBeam>::set_beam)
      .add_property("detector", &Experiment<MonochromaticBeam>::get_detector, &Experiment<MonochromaticBeam>::set_detector)
      .add_property(
        "goniometer", &Experiment<MonochromaticBeam>::get_goniometer, &Experiment<MonochromaticBeam>::set_goniometer)
      .add_property("scan", &Experiment<MonochromaticBeam>::get_scan, &Experiment<MonochromaticBeam>::set_scan)
      .add_property("crystal", &Experiment<MonochromaticBeam>::get_crystal, &Experiment<MonochromaticBeam>::set_crystal)
      .add_property("profile", &Experiment<MonochromaticBeam>::get_profile, &Experiment<MonochromaticBeam>::set_profile)
      .add_property("imageset", &Experiment<MonochromaticBeam>::get_imageset, &Experiment<MonochromaticBeam>::set_imageset)
      .add_property(
        "scaling_model", &Experiment<MonochromaticBeam>::get_scaling_model, &Experiment<MonochromaticBeam>::set_scaling_model)
      .add_property(
        "identifier", &Experiment<MonochromaticBeam>::get_identifier, &Experiment<MonochromaticBeam>::set_identifier)
      .def("__contains__", experiment_contains_pointers<MonochromaticBeam>::beam())
      .def("__contains__", experiment_contains_pointers<MonochromaticBeam>::detector())
      .def("__contains__", experiment_contains_pointers<MonochromaticBeam>::goniometer())
      .def("__contains__", experiment_contains_pointers<MonochromaticBeam>::scan())
      .def("__contains__", experiment_contains_pointers<MonochromaticBeam>::crystal())
      .def("__contains__", experiment_contains_pointers<MonochromaticBeam>::object())
      .def("__eq__", &Experiment<MonochromaticBeam>::operator==)
      .def("is_consistent", &Experiment<MonochromaticBeam>::is_consistent)
      .def("is_still",
           &Experiment<MonochromaticBeam>::is_still,
           "Check if this experiment represents a still image")
      .def("is_sequence",
           &Experiment<MonochromaticBeam>::is_sequence,
           "Check if this experiment represents swept rotation image(s)")
      .def_pickle(ExperimentPickleSuite<MonochromaticBeam>());
    class_<Experiment<TOFBeam> >("Experiment")
      .def(init<boost::shared_ptr<TOFBeam>,
                boost::shared_ptr<Detector>,
                boost::shared_ptr<Goniometer>,
                boost::shared_ptr<Scan>,
                boost::shared_ptr<CrystalBase>,
                boost::python::object,
                boost::python::object,
                boost::python::object,
                std::string>((arg("beam") = boost::shared_ptr<TOFBeam>(),
                              arg("detector") = boost::shared_ptr<Detector>(),
                              arg("goniometer") = boost::shared_ptr<Goniometer>(),
                              arg("scan") = boost::shared_ptr<Scan>(),
                              arg("crystal") = boost::shared_ptr<CrystalBase>(),
                              arg("profile") = boost::python::object(),
                              arg("imageset") = boost::python::object(),
                              arg("scaling_model") = boost::python::object(),
                              arg("identifier") = "")))
      .add_property("beam", &Experiment<TOFBeam>::get_beam, &Experiment<TOFBeam>::set_beam)
      .add_property("detector", &Experiment<TOFBeam>::get_detector, &Experiment<TOFBeam>::set_detector)
      .add_property(
        "goniometer", &Experiment<TOFBeam>::get_goniometer, &Experiment<TOFBeam>::set_goniometer)
      .add_property("scan", &Experiment<TOFBeam>::get_scan, &Experiment<TOFBeam>::set_scan)
      .add_property("crystal", &Experiment<TOFBeam>::get_crystal, &Experiment<TOFBeam>::set_crystal)
      .add_property("profile", &Experiment<TOFBeam>::get_profile, &Experiment<TOFBeam>::set_profile)
      .add_property("imageset", &Experiment<TOFBeam>::get_imageset, &Experiment<TOFBeam>::set_imageset)
      .add_property(
        "scaling_model", &Experiment<TOFBeam>::get_scaling_model, &Experiment<TOFBeam>::set_scaling_model)
      .add_property(
        "identifier", &Experiment<TOFBeam>::get_identifier, &Experiment<TOFBeam>::set_identifier)
      .def("__contains__", experiment_contains_pointers<TOFBeam>::beam())
      .def("__contains__", experiment_contains_pointers<TOFBeam>::detector())
      .def("__contains__", experiment_contains_pointers<TOFBeam>::goniometer())
      .def("__contains__", experiment_contains_pointers<TOFBeam>::scan())
      .def("__contains__", experiment_contains_pointers<TOFBeam>::crystal())
      .def("__contains__", experiment_contains_pointers<TOFBeam>::object())
      .def("__eq__", &Experiment<TOFBeam>::operator==)
      .def("is_consistent", &Experiment<TOFBeam>::is_consistent)
      .def("is_still",
           &Experiment<TOFBeam>::is_still,
           "Check if this experiment represents a still image")
      .def("is_sequence",
           &Experiment<TOFBeam>::is_sequence,
           "Check if this experiment represents swept rotation image(s)")
      .def_pickle(ExperimentPickleSuite<TOFBeam>());
  }

}}}  // namespace dxtbx::model::boost_python
