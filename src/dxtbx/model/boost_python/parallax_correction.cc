/*
 * parallax_correction.cc
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
#include <dxtbx/model/parallax_correction.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  // Parallax functions are overloaded so declare the forms for boost
  typedef vec2<double> (*parallax_orig_type)(double d,
                                             double la,
                                             vec2<double> xy0,
                                             vec2<double> xy);
  typedef vec2<double> (*parallax_axis_type)(double mu,
                                             double t0,
                                             vec2<double> xy,
                                             vec3<double> fast_axis,
                                             vec3<double> slow_axis,
                                             vec3<double> origin);

  void export_parallax_correction() {
    def("parallax_correction",
        (parallax_orig_type)&parallax_correction,
        (arg("d"), arg("la"), arg("xy0"), arg("xy")));
    def("parallax_correction",
        (parallax_axis_type)&parallax_correction,
        (arg("mu"), arg("t0"), arg("xy"), arg("fast"), arg("slow"), arg("origin")));
    def("parallax_correction_inv",
        (parallax_orig_type)&parallax_correction_inv,
        (arg("d"), arg("la"), arg("xy0"), arg("xy")));
    def("parallax_correction_inv",
        (parallax_axis_type)&parallax_correction_inv,
        (arg("mu"), arg("t0"), arg("xy"), arg("fast"), arg("slow"), arg("origin")));
  }

}}}  // namespace dxtbx::model::boost_python
