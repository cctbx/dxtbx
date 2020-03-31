/*
 * crystal.cc
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
#include <boost_adaptbx/optional_conversions.h>
#include <string>
#include <sstream>
#include <dxtbx/model/crystal.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using boost_adaptbx::optional_conversions::to_and_from_python;

  static Crystal *make_crystal_default(const vec3<double> &real_space_a,
                                       const vec3<double> &real_space_b,
                                       const vec3<double> &real_space_c,
                                       const cctbx::sgtbx::space_group &space_group) {
    Crystal *crystal =
      new Crystal(real_space_a, real_space_b, real_space_c, space_group);

    return crystal;
  }

  static Crystal *make_crystal_with_symbol(const vec3<double> &real_space_a,
                                           const vec3<double> &real_space_b,
                                           const vec3<double> &real_space_c,
                                           const std::string &space_group_symbol) {
    Crystal *crystal = new Crystal(
      real_space_a,
      real_space_b,
      real_space_c,
      cctbx::sgtbx::space_group(cctbx::sgtbx::space_group_symbols(space_group_symbol)));

    return crystal;
  }

  static Crystal *make_crystal_with_A(const mat3<double> &A,
                                      const cctbx::sgtbx::space_group &space_group,
                                      const bool &reciprocal) {
    Crystal *crystal = new Crystal(A, space_group, reciprocal);

    return crystal;
  }

  static Crystal *make_crystal_with_A_symbol(const mat3<double> &A,
                                             const std::string &space_group_symbol,
                                             const bool &reciprocal) {
    Crystal *crystal = new Crystal(
      A,
      cctbx::sgtbx::space_group(cctbx::sgtbx::space_group_symbols(space_group_symbol)),
      reciprocal);

    return crystal;
  }

  static MosaicCrystalKabsch2010 *make_kabsch2010_mosaic_crystal_default(
    const vec3<double> &real_space_a,
    const vec3<double> &real_space_b,
    const vec3<double> &real_space_c,
    const cctbx::sgtbx::space_group &space_group) {
    MosaicCrystalKabsch2010 *crystal = new MosaicCrystalKabsch2010(
      real_space_a, real_space_b, real_space_c, space_group);

    return crystal;
  }

  static MosaicCrystalKabsch2010 *make_kabsch2010_mosaic_crystal_with_symbol(
    const vec3<double> &real_space_a,
    const vec3<double> &real_space_b,
    const vec3<double> &real_space_c,
    const std::string &space_group_symbol) {
    MosaicCrystalKabsch2010 *crystal = new MosaicCrystalKabsch2010(
      real_space_a,
      real_space_b,
      real_space_c,
      cctbx::sgtbx::space_group(cctbx::sgtbx::space_group_symbols(space_group_symbol)));

    return crystal;
  }

  static MosaicCrystalSauter2014 *make_sauter2014_mosaic_crystal_default(
    const vec3<double> &real_space_a,
    const vec3<double> &real_space_b,
    const vec3<double> &real_space_c,
    const cctbx::sgtbx::space_group &space_group) {
    MosaicCrystalSauter2014 *crystal = new MosaicCrystalSauter2014(
      real_space_a, real_space_b, real_space_c, space_group);

    return crystal;
  }

  static MosaicCrystalSauter2014 *make_sauter2014_mosaic_crystal_with_symbol(
    const vec3<double> &real_space_a,
    const vec3<double> &real_space_b,
    const vec3<double> &real_space_c,
    const std::string &space_group_symbol) {
    MosaicCrystalSauter2014 *crystal = new MosaicCrystalSauter2014(
      real_space_a,
      real_space_b,
      real_space_c,
      cctbx::sgtbx::space_group(cctbx::sgtbx::space_group_symbols(space_group_symbol)));

    return crystal;
  }

  static void Crystal_set_A_at_scan_points_from_tuple(CrystalBase &self,
                                                      boost::python::tuple l) {
    scitbx::af::shared<mat3<double> > A_list;
    for (std::size_t i = 0; i < boost::python::len(l); ++i) {
      mat3<double> A = boost::python::extract<mat3<double> >(l[i]);
      A_list.push_back(A);
    }
    self.set_A_at_scan_points(A_list.const_ref());
  }

  static void Crystal_set_A_at_scan_points_from_list(CrystalBase &self,
                                                     boost::python::list l) {
    scitbx::af::shared<mat3<double> > A_list;
    for (std::size_t i = 0; i < boost::python::len(l); ++i) {
      mat3<double> A = boost::python::extract<mat3<double> >(l[i]);
      A_list.push_back(A);
    }
    self.set_A_at_scan_points(A_list.const_ref());
  }

  static void Crystal_set_B_covariance_from_tuple(CrystalBase &self,
                                                  boost::python::object obj) {
    scitbx::af::versa<double, scitbx::af::c_grid<2> > B_cov(
      scitbx::af::c_grid<2>(9, 9));
    DXTBX_ASSERT(boost::python::len(obj) == 9 * 9);
    for (std::size_t i = 0; i < boost::python::len(obj); ++i) {
      B_cov[i] = boost::python::extract<double>(obj[i]);
    }
    self.set_B_covariance(B_cov.const_ref());
  }

  struct CrystalPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(const Crystal &obj) {
      scitbx::af::shared<vec3<double> > real_space_v = obj.get_real_space_vectors();
      return boost::python::make_tuple(
        real_space_v[0], real_space_v[1], real_space_v[2], obj.get_space_group());
    }

    static boost::python::tuple getstate(boost::python::object obj) {
      const Crystal &crystal = boost::python::extract<const Crystal &>(obj)();
      return boost::python::make_tuple(obj.attr("__dict__"),
                                       crystal.get_A_at_scan_points(),
                                       crystal.get_B_covariance(),
                                       crystal.get_B_covariance_at_scan_points(),
                                       crystal.get_recalculated_unit_cell(),
                                       crystal.get_recalculated_cell_parameter_sd(),
                                       crystal.get_recalculated_cell_volume_sd());
    }

    static void setstate(boost::python::object obj, boost::python::tuple state) {
      Crystal &crystal = boost::python::extract<Crystal &>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 7);

      // restore the object's __dict__
      boost::python::dict d =
        boost::python::extract<boost::python::dict>(obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      scitbx::af::const_ref<mat3<double> > A_scanpoints =
        boost::python::extract<scitbx::af::const_ref<mat3<double> > >(state[1]);
      scitbx::af::const_ref<double, scitbx::af::c_grid<2> > cov_B =
        boost::python::extract<scitbx::af::const_ref<double, scitbx::af::c_grid<2> > >(
          state[2]);
      scitbx::af::const_ref<double, scitbx::af::c_grid<3> > cov_B_scanpoints =
        boost::python::extract<scitbx::af::const_ref<double, scitbx::af::c_grid<3> > >(
          state[3]);
      boost::optional<cctbx::uctbx::unit_cell> recalculated_unit_cell =
        boost::python::extract<boost::optional<cctbx::uctbx::unit_cell> >(state[4]);
      scitbx::af::small<double, 6> recalculated_cell_parameter_sd =
        boost::python::extract<scitbx::af::small<double, 6> >(state[5]);
      double recalculated_cell_volume_sd = boost::python::extract<double>(state[6]);
      crystal.set_A_at_scan_points(A_scanpoints);
      crystal.set_B_covariance(cov_B);
      crystal.set_B_covariance_at_scan_points(cov_B_scanpoints);
      if (recalculated_unit_cell) {
        crystal.set_recalculated_unit_cell(*recalculated_unit_cell);
      }
      crystal.set_recalculated_cell_parameter_sd(recalculated_cell_parameter_sd);
      crystal.set_recalculated_cell_volume_sd(recalculated_cell_volume_sd);
    }

    static bool getstate_manages_dict() {
      return true;
    }
  };

  struct MosaicCrystalKabsch2010PickleSuite : CrystalPickleSuite {
    static boost::python::tuple getinitargs(const MosaicCrystalKabsch2010 &obj) {
      scitbx::af::shared<vec3<double> > real_space_v = obj.get_real_space_vectors();
      return boost::python::make_tuple(
        real_space_v[0], real_space_v[1], real_space_v[2], obj.get_space_group());
    }

    static boost::python::tuple getstate(boost::python::object obj) {
      const MosaicCrystalKabsch2010 &crystal =
        boost::python::extract<const MosaicCrystalKabsch2010 &>(obj)();
      return boost::python::make_tuple(obj.attr("__dict__"),
                                       crystal.get_A_at_scan_points(),
                                       crystal.get_B_covariance(),
                                       crystal.get_B_covariance_at_scan_points(),
                                       crystal.get_recalculated_unit_cell(),
                                       crystal.get_recalculated_cell_parameter_sd(),
                                       crystal.get_recalculated_cell_volume_sd(),
                                       crystal.get_mosaicity());
    }

    static void setstate(boost::python::object obj, boost::python::tuple state) {
      MosaicCrystalKabsch2010 &crystal =
        boost::python::extract<MosaicCrystalKabsch2010 &>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 8);

      // restore the object's __dict__
      boost::python::dict d =
        boost::python::extract<boost::python::dict>(obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      scitbx::af::const_ref<mat3<double> > A_scanpoints =
        boost::python::extract<scitbx::af::const_ref<mat3<double> > >(state[1]);
      scitbx::af::const_ref<double, scitbx::af::c_grid<2> > cov_B =
        boost::python::extract<scitbx::af::const_ref<double, scitbx::af::c_grid<2> > >(
          state[2]);
      scitbx::af::const_ref<double, scitbx::af::c_grid<3> > cov_B_scanpoints =
        boost::python::extract<scitbx::af::const_ref<double, scitbx::af::c_grid<3> > >(
          state[3]);
      boost::optional<cctbx::uctbx::unit_cell> recalculated_unit_cell =
        boost::python::extract<boost::optional<cctbx::uctbx::unit_cell> >(state[4]);
      scitbx::af::small<double, 6> recalculated_cell_parameter_sd =
        boost::python::extract<scitbx::af::small<double, 6> >(state[5]);
      double recalculated_cell_volume_sd = boost::python::extract<double>(state[6]);
      double mosaicity = boost::python::extract<double>(state[7]);
      crystal.set_A_at_scan_points(A_scanpoints);
      crystal.set_B_covariance(cov_B);
      crystal.set_B_covariance_at_scan_points(cov_B_scanpoints);
      if (recalculated_unit_cell) {
        crystal.set_recalculated_unit_cell(*recalculated_unit_cell);
      }
      crystal.set_recalculated_cell_parameter_sd(recalculated_cell_parameter_sd);
      crystal.set_recalculated_cell_volume_sd(recalculated_cell_volume_sd);
      crystal.set_mosaicity(mosaicity);
    }
  };

  struct MosaicCrystalSauter2014PickleSuite : CrystalPickleSuite {
    static boost::python::tuple getinitargs(const MosaicCrystalSauter2014 &obj) {
      scitbx::af::shared<vec3<double> > real_space_v = obj.get_real_space_vectors();
      return boost::python::make_tuple(
        real_space_v[0], real_space_v[1], real_space_v[2], obj.get_space_group());
    }

    static boost::python::tuple getstate(boost::python::object obj) {
      const MosaicCrystalSauter2014 &crystal =
        boost::python::extract<const MosaicCrystalSauter2014 &>(obj)();
      return boost::python::make_tuple(obj.attr("__dict__"),
                                       crystal.get_A_at_scan_points(),
                                       crystal.get_B_covariance(),
                                       crystal.get_B_covariance_at_scan_points(),
                                       crystal.get_recalculated_unit_cell(),
                                       crystal.get_recalculated_cell_parameter_sd(),
                                       crystal.get_recalculated_cell_volume_sd(),
                                       crystal.get_half_mosaicity_deg(),
                                       crystal.get_domain_size_ang());
    }

    static void setstate(boost::python::object obj, boost::python::tuple state) {
      MosaicCrystalSauter2014 &crystal =
        boost::python::extract<MosaicCrystalSauter2014 &>(obj)();
      DXTBX_ASSERT(boost::python::len(state) == 9);

      // restore the object's __dict__
      boost::python::dict d =
        boost::python::extract<boost::python::dict>(obj.attr("__dict__"))();
      d.update(state[0]);

      // restore the internal state of the C++ object
      scitbx::af::const_ref<mat3<double> > A_scanpoints =
        boost::python::extract<scitbx::af::const_ref<mat3<double> > >(state[1]);
      scitbx::af::const_ref<double, scitbx::af::c_grid<2> > cov_B =
        boost::python::extract<scitbx::af::const_ref<double, scitbx::af::c_grid<2> > >(
          state[2]);
      scitbx::af::const_ref<double, scitbx::af::c_grid<3> > cov_B_scanpoints =
        boost::python::extract<scitbx::af::const_ref<double, scitbx::af::c_grid<3> > >(
          state[3]);
      boost::optional<cctbx::uctbx::unit_cell> recalculated_unit_cell =
        boost::python::extract<boost::optional<cctbx::uctbx::unit_cell> >(state[4]);
      scitbx::af::small<double, 6> recalculated_cell_parameter_sd =
        boost::python::extract<scitbx::af::small<double, 6> >(state[5]);
      double recalculated_cell_volume_sd = boost::python::extract<double>(state[6]);
      double half_mosaicity_deg = boost::python::extract<double>(state[7]);
      double domain_size_ang = boost::python::extract<double>(state[8]);
      crystal.set_A_at_scan_points(A_scanpoints);
      crystal.set_B_covariance(cov_B);
      crystal.set_B_covariance_at_scan_points(cov_B_scanpoints);
      if (recalculated_unit_cell) {
        crystal.set_recalculated_unit_cell(*recalculated_unit_cell);
      }
      crystal.set_recalculated_cell_parameter_sd(recalculated_cell_parameter_sd);
      crystal.set_recalculated_cell_volume_sd(recalculated_cell_volume_sd);
      crystal.set_half_mosaicity_deg(half_mosaicity_deg);
      crystal.set_domain_size_ang(domain_size_ang);
    }
  };

  inline void CrystalBase_set_unit_cell_real_space_vectors(
    CrystalBase &self,
    scitbx::vec3<double> &real_space_a,
    scitbx::vec3<double> &real_space_b,
    scitbx::vec3<double> &real_space_c) {
    self.set_unit_cell(real_space_a, real_space_b, real_space_c);
  }

  inline void CrystalBase_set_unit_cell(CrystalBase &self,
                                        cctbx::uctbx::unit_cell &unit_cell) {
    self.set_unit_cell(unit_cell);
  }

  inline void CrystalBase_set_recalculated_unit_cell(
    CrystalBase &self,
    cctbx::uctbx::unit_cell &unit_cell) {
    self.set_recalculated_unit_cell(unit_cell);
  }

  inline void CrystalBase_set_recalculated_cell_parameter_sd(
    CrystalBase &self,
    const scitbx::af::small<double, 6> &unit_cell_sd) {
    self.set_recalculated_cell_parameter_sd(unit_cell_sd);
  }

  void export_crystal() {
    // Expose the optional values
    to_and_from_python<boost::optional<cctbx::uctbx::unit_cell> >();

    class_<CrystalBase, boost::noncopyable>("CrystalBase", no_init)
      .def("set_unit_cell", CrystalBase_set_unit_cell_real_space_vectors)
      .def("set_unit_cell", CrystalBase_set_unit_cell)
      .def("update_B", &CrystalBase::update_B)
      .def("set_U", &CrystalBase::set_U)
      .def("get_U", &CrystalBase::get_U)
      .def("set_B", &CrystalBase::set_B)
      .def("get_B", &CrystalBase::get_B)
      .def("set_A", &CrystalBase::set_A)
      .def("get_A", &CrystalBase::get_A)
      .def("get_unit_cell", &CrystalBase::get_unit_cell)
      .def("get_real_space_vectors", &CrystalBase::get_real_space_vectors)
      .def("set_space_group", &CrystalBase::set_space_group)
      .def("get_space_group", &CrystalBase::get_space_group)
      .add_property("num_scan_points", &CrystalBase::get_num_scan_points)
      .def("get_num_scan_points", &CrystalBase::get_num_scan_points)
      .def("set_A_at_scan_points", &CrystalBase::set_A_at_scan_points)
      .def("set_A_at_scan_points", &Crystal_set_A_at_scan_points_from_tuple)
      .def("set_A_at_scan_points", &Crystal_set_A_at_scan_points_from_list)
      .def("get_A_at_scan_point", &CrystalBase::get_A_at_scan_point)
      .def("get_B_at_scan_point", &CrystalBase::get_B_at_scan_point)
      .def("get_U_at_scan_point", &CrystalBase::get_U_at_scan_point)
      .def("get_unit_cell_at_scan_point", &CrystalBase::get_unit_cell_at_scan_point)
      .def("reset_scan_points", &CrystalBase::reset_scan_points)
      .def("change_basis", &CrystalBase::change_basis)
      .def("update", &CrystalBase::update)
      .def("rotate_around_origin",
           &CrystalBase::rotate_around_origin,
           (arg("axis"), arg("angle"), arg("deg") = true))
      .def("is_similar_to",
           &CrystalBase::is_similar_to,
           (arg("other"),
            arg("angle_tolerance") = 0.01,
            arg("uc_rel_length_tolerance") = 0.01,
            arg("uc_abs_angle_tolerance") = 1.0))
      .def("get_B_covariance", &CrystalBase::get_B_covariance)
      .def("set_B_covariance", &CrystalBase::set_B_covariance)
      .def("set_B_covariance", &Crystal_set_B_covariance_from_tuple)
      .def("set_B_covariance_at_scan_points",
           &CrystalBase::set_B_covariance_at_scan_points)
      .def("get_B_covariance_at_scan_point",
           &CrystalBase::get_B_covariance_at_scan_point)
      .def("get_B_covariance_at_scan_points",
           &CrystalBase::get_B_covariance_at_scan_points)
      .def("get_cell_parameter_sd", &CrystalBase::get_cell_parameter_sd)
      .def("get_cell_volume_sd", &CrystalBase::get_cell_volume_sd)
      .def("get_recalculated_cell_volume_sd",
           &CrystalBase::get_recalculated_cell_volume_sd)
      .def("set_recalculated_cell_volume_sd",
           &CrystalBase::set_recalculated_cell_volume_sd)
      .def("get_cell_parameter_sd_at_scan_point",
           &CrystalBase::get_cell_parameter_sd_at_scan_point)
      .def("reset_unit_cell_errors", &CrystalBase::reset_unit_cell_errors)
      .def("set_recalculated_unit_cell", &CrystalBase_set_recalculated_unit_cell)
      .def("get_recalculated_unit_cell", &CrystalBase::get_recalculated_unit_cell)
      .def("set_recalculated_cell_parameter_sd",
           &CrystalBase_set_recalculated_cell_parameter_sd)
      .def("get_recalculated_cell_parameter_sd",
           &CrystalBase::get_recalculated_cell_parameter_sd)
      .def("__eq__", &CrystalBase::operator==)
      .def("__ne__", &CrystalBase::operator!=);

    class_<Crystal, bases<CrystalBase> >("Crystal", no_init)
      .def(init<const Crystal &>())
      .def("__init__",
           make_constructor(&make_crystal_default,
                            default_call_policies(),
                            (arg("real_space_a"),
                             arg("real_space_b"),
                             arg("real_space_c"),
                             arg("space_group"))))
      .def("__init__",
           make_constructor(&make_crystal_with_symbol,
                            default_call_policies(),
                            (arg("real_space_a"),
                             arg("real_space_b"),
                             arg("real_space_c"),
                             arg("space_group_symbol"))))
      .def("__init__",
           make_constructor(&make_crystal_with_A,
                            default_call_policies(),
                            (arg("A"), arg("space_group"), arg("reciprocal") = true)))
      .def("__init__",
           make_constructor(
             &make_crystal_with_A_symbol,
             default_call_policies(),
             (arg("A"), arg("space_group_symbol"), arg("reciprocal") = true)))
      .def_pickle(CrystalPickleSuite());

    // Create member-function pointers to specific is_similar_to overloads
    // - each of these crystals has a custom implementation along with the
    //   inherited interface, and we want to expose these explicitly
    bool (MosaicCrystalKabsch2010::*kabsch_is_similar_to)(
      const CrystalBase &, double, double, double, double) const =
      &MosaicCrystalKabsch2010::is_similar_to;
    bool (MosaicCrystalSauter2014::*sauter_is_similar_to)(
      const CrystalBase &, double, double, double, double, double) const =
      &MosaicCrystalSauter2014::is_similar_to;

    class_<MosaicCrystalKabsch2010, bases<Crystal> >("MosaicCrystalKabsch2010", no_init)
      .def(init<const MosaicCrystalKabsch2010 &>())
      .def(init<const Crystal &>())
      .def("__init__",
           make_constructor(&make_kabsch2010_mosaic_crystal_default,
                            default_call_policies(),
                            (arg("real_space_a"),
                             arg("real_space_b"),
                             arg("real_space_c"),
                             arg("space_group"))))
      .def("__init__",
           make_constructor(&make_kabsch2010_mosaic_crystal_with_symbol,
                            default_call_policies(),
                            (arg("real_space_a"),
                             arg("real_space_b"),
                             arg("real_space_c"),
                             arg("space_group_symbol"))))
      .def("is_similar_to",
           kabsch_is_similar_to,
           (arg("other"),
            arg("angle_tolerance") = 0.01,
            arg("uc_rel_length_tolerance") = 0.01,
            arg("uc_abs_angle_tolerance") = 1.0,
            arg("mosaicity_tolerance") = 0.8))
      .def(
        "get_mosaicity", &MosaicCrystalKabsch2010::get_mosaicity, (arg("deg") = true))
      .def("set_mosaicity",
           &MosaicCrystalKabsch2010::set_mosaicity,
           (arg("mosaicity"), arg("deg") = true))
      .def_pickle(MosaicCrystalKabsch2010PickleSuite());

    class_<MosaicCrystalSauter2014, bases<Crystal> >("MosaicCrystalSauter2014", no_init)
      .def(init<const MosaicCrystalSauter2014 &>())
      .def(init<const Crystal &>())
      .def("__init__",
           make_constructor(&make_sauter2014_mosaic_crystal_default,
                            default_call_policies(),
                            (arg("real_space_a"),
                             arg("real_space_b"),
                             arg("real_space_c"),
                             arg("space_group"))))
      .def("__init__",
           make_constructor(&make_sauter2014_mosaic_crystal_with_symbol,
                            default_call_policies(),
                            (arg("real_space_a"),
                             arg("real_space_b"),
                             arg("real_space_c"),
                             arg("space_group_symbol"))))
      .def("is_similar_to",
           sauter_is_similar_to,
           (arg("other"),
            arg("angle_tolerance") = 0.01,
            arg("uc_rel_length_tolerance") = 0.01,
            arg("uc_abs_angle_tolerance") = 1.0,
            arg("half_mosaicity_tolerance") = 0.4,
            arg("domain_size_tolerance") = 1.0))
      .def("get_half_mosaicity_deg", &MosaicCrystalSauter2014::get_half_mosaicity_deg)
      .def("set_half_mosaicity_deg",
           &MosaicCrystalSauter2014::set_half_mosaicity_deg,
           (arg("half_mosaicity_deg")))
      .def("get_domain_size_ang", &MosaicCrystalSauter2014::get_domain_size_ang)
      .def("set_domain_size_ang",
           &MosaicCrystalSauter2014::set_domain_size_ang,
           (arg("domain_size_ang")))
      .def_pickle(MosaicCrystalSauter2014PickleSuite());

    register_ptr_to_python<boost::shared_ptr<CrystalBase> >();
    register_ptr_to_python<boost::shared_ptr<Crystal> >();
    register_ptr_to_python<boost::shared_ptr<MosaicCrystalKabsch2010> >();
    register_ptr_to_python<boost::shared_ptr<MosaicCrystalSauter2014> >();
  }

}}}  // namespace dxtbx::model::boost_python
