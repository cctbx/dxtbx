/*
 * goniometer.cc
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
#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <sstream>
#include <scitbx/constants.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/kappa_goniometer.h>
#include <dxtbx/model/multi_axis_goniometer.h>
#include <dxtbx/model/boost_python/to_from_dict.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;

  static void rotate_around_origin(Goniometer &goniometer,
                                   vec3<double> axis,
                                   double angle,
                                   bool deg) {
    double angle_rad = deg ? deg_as_rad(angle) : angle;
    goniometer.rotate_around_origin(axis, angle_rad);
  }

  std::string goniometer_to_string(const Goniometer &goniometer) {
    std::stringstream ss;
    ss << goniometer;
    return ss.str();
  }

  struct GoniometerPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(const Goniometer &obj) {
      return boost::python::make_tuple(obj.get_rotation_axis_datum(),
                                       obj.get_fixed_rotation(),
                                       obj.get_setting_rotation());
    }

    static boost::python::tuple getstate(boost::python::object obj) {
      const Goniometer &goniometer = boost::python::extract<const Goniometer &>(obj)();
      return boost::python::make_tuple(
        obj.attr("__dict__"), goniometer.get_setting_rotation_at_scan_points());
    }

    static void setstate(boost::python::object obj, boost::python::tuple state) {
      Goniometer &goniometer = boost::python::extract<Goniometer &>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 2);

      // restore the object's __dict__
      boost::python::dict d =
        boost::python::extract<boost::python::dict>(obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      scitbx::af::const_ref<mat3<double> > S_list =
        boost::python::extract<scitbx::af::const_ref<mat3<double> > >(state[1]);
      goniometer.set_setting_rotation_at_scan_points(S_list);
    }

    static bool getstate_manages_dict() {
      return true;
    }
  };

  static void Goniometer_set_S_at_scan_points_from_tuple(Goniometer &goniometer,
                                                         boost::python::tuple l) {
    scitbx::af::shared<mat3<double> > S_list;
    for (std::size_t i = 0; i < boost::python::len(l); ++i) {
      mat3<double> S = boost::python::extract<mat3<double> >(l[i]);
      S_list.push_back(S);
    }
    goniometer.set_setting_rotation_at_scan_points(S_list.const_ref());
  }

  static void Goniometer_set_S_at_scan_points_from_list(Goniometer &goniometer,
                                                        boost::python::list l) {
    scitbx::af::shared<mat3<double> > S_list;
    for (std::size_t i = 0; i < boost::python::len(l); ++i) {
      mat3<double> S = boost::python::extract<mat3<double> >(l[i]);
      S_list.push_back(S);
    }
    goniometer.set_setting_rotation_at_scan_points(S_list.const_ref());
  }

  template <>
  boost::python::dict to_dict<Goniometer>(const Goniometer &obj) {
    boost::python::dict result;
    result["rotation_axis"] = obj.get_rotation_axis_datum();
    result["fixed_rotation"] = obj.get_fixed_rotation();
    result["setting_rotation"] = obj.get_setting_rotation();
    if (obj.get_num_scan_points() > 0) {
      boost::python::list l;
      scitbx::af::shared<mat3<double> > setting_rotation_at_scan_points =
        obj.get_setting_rotation_at_scan_points();
      for (scitbx::af::shared<mat3<double> >::iterator it =
             setting_rotation_at_scan_points.begin();
           it != setting_rotation_at_scan_points.end();
           ++it) {
        l.append(boost::python::tuple(*it));
      }
      result["setting_rotation_at_scan_points"] = l;
    }
    return result;
  }

  template <>
  Goniometer *from_dict<Goniometer>(boost::python::dict obj) {
    Goniometer *g = new Goniometer(
      boost::python::extract<vec3<double> >(obj["rotation_axis"]),
      boost::python::extract<mat3<double> >(obj.get(
        "fixed_rotation", mat3<double>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))),
      boost::python::extract<mat3<double> >(
        obj.get("setting_rotation",
                mat3<double>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))));
    if (obj.has_key("setting_rotation_at_scan_points")) {
      boost::python::list S_at_scan_points =
        boost::python::extract<boost::python::list>(
          obj["setting_rotation_at_scan_points"]);
      Goniometer_set_S_at_scan_points_from_list(*g, S_at_scan_points);
    }
    return g;
  }

  void export_goniometer() {
    class_<GoniometerBase>("GoniometerBase");

    class_<Goniometer, boost::shared_ptr<Goniometer>, bases<GoniometerBase> >(
      "Goniometer")
      .def(init<const Goniometer &>())
      .def(init<vec3<double> >((arg("rotation_axis"))))
      .def(init<vec3<double>, mat3<double> >(
        (arg("rotation_axis"), arg("fixed_rotation_matrix"))))
      .def(init<vec3<double>, mat3<double>, mat3<double> >(
        (arg("rotation_axis"),
         arg("fixed_rotation_matrix"),
         arg("setting_rotation_matrix"))))
      .def("get_rotation_axis", &Goniometer::get_rotation_axis)
      .def("set_rotation_axis", &Goniometer::set_rotation_axis)
      .def("get_rotation_axis_datum", &Goniometer::get_rotation_axis_datum)
      .def("set_rotation_axis_datum", &Goniometer::set_rotation_axis_datum)
      .def("get_fixed_rotation", &Goniometer::get_fixed_rotation)
      .def("set_fixed_rotation", &Goniometer::set_fixed_rotation)
      .def("get_setting_rotation", &Goniometer::get_setting_rotation)
      .def("set_setting_rotation", &Goniometer::set_setting_rotation)
      .add_property("num_scan_points", &Goniometer::get_num_scan_points)
      .def("get_num_scan_points", &Goniometer::get_num_scan_points)
      .def("set_setting_rotation_at_scan_points",
           &Goniometer::set_setting_rotation_at_scan_points)
      .def("set_setting_rotation_at_scan_points",
           &Goniometer_set_S_at_scan_points_from_tuple)
      .def("set_setting_rotation_at_scan_points",
           &Goniometer_set_S_at_scan_points_from_list)
      .def("get_setting_rotation_at_scan_points",
           &Goniometer::get_setting_rotation_at_scan_points)
      .def("get_setting_rotation_at_scan_point",
           &Goniometer::get_setting_rotation_at_scan_point)
      .def("reset_scan_points", &Goniometer::reset_scan_points)
      .def("rotate_around_origin",
           &rotate_around_origin,
           (arg("axis"), arg("angle"), arg("deg") = true))
      .def("__eq__", &Goniometer::operator==)
      .def("__ne__", &Goniometer::operator!=)
      .def("is_similar_to",
           &Goniometer::is_similar_to,
           (arg("other"),
            arg("rotation_axis_tolerance") = 1e-6,
            arg("fixed_rotation_tolerance") = 1e-6,
            arg("setting_rotation_tolerance") = 1e-6))
      .def("__str__", &goniometer_to_string)
      .def("to_dict", &to_dict<Goniometer>)
      .def(
        "from_dict", &from_dict<Goniometer>, return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(GoniometerPickleSuite());
  }

}}}  // namespace dxtbx::model::boost_python
