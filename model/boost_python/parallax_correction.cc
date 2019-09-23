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

  vec2<double> parallax_correction_original(double d,
                                            double la,
                                            vec2<double> xy0,
                                            vec2<double> xy) {
    return parallax_correction(d, la, xy0, xy);
  }

  vec2<double> parallax_correction_with_panel_data(double mu,
                                                   double t0,
                                                   vec2<double> xy,
                                                   vec3<double> fast,
                                                   vec3<double> slow,
                                                   vec3<double> origin) {
    return parallax_correction(mu, t0, xy, fast, slow, origin);
  }

  void export_parallax_correction() {
    def("parallax_correction",
        &parallax_correction_original,
        (arg("d"), arg("la"), arg("xy0"), arg("xy")));
    def("parallax_correction",
        &parallax_correction_with_panel_data,
        (arg("mu"), arg("t0"), arg("xy"), arg("fast"), arg("slow"), arg("origin")));
    def("parallax_correction_inv",
        &parallax_correction_inv,
        (arg("d"), arg("la"), arg("xy0"), arg("xy")));
    def("parallax_correction_inv2",
        &parallax_correction_inv2,
        (arg("mu"), arg("t0"), arg("xy"), arg("fast"), arg("slow"), arg("origin")));
  }

}}}  // namespace dxtbx::model::boost_python
