/*
 * experiment_ext.cc
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

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_beam();
  void export_goniometer();
  void export_kappa_goniometer();
  void export_multi_axis_goniometer();
  void export_panel();
  void export_detector();
  void export_scan();
  void export_scan_helpers();
  void export_crystal();
  void export_parallax_correction();
  void export_pixel_to_millimeter();
  void export_experiment();
  void export_experiment_list();
  void export_spectrum();

  BOOST_PYTHON_MODULE(dxtbx_model_ext) {
    export_beam();
    export_goniometer();
    export_kappa_goniometer();
    export_multi_axis_goniometer();
    export_panel();
    export_detector();
    export_scan();
    export_scan_helpers();
    export_crystal();
    export_parallax_correction();
    export_pixel_to_millimeter();
    export_experiment();
    export_experiment_list();
    export_spectrum();
  }

}}}  // namespace dxtbx::model::boost_python
