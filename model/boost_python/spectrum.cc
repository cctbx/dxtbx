/*
 * spectrum.cc
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <string>
#include <sstream>
#include <dxtbx/model/spectrum.h>
#include <dxtbx/model/boost_python/to_from_dict.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;
  using scitbx::rad_as_deg;

  std::string spectrum_to_string(const Spectrum &spectrum) {
    std::stringstream ss;
    ss << spectrum;
    return ss.str();
  }

  struct SpectrumPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(const Spectrum &obj) {
      return boost::python::make_tuple(obj.get_energies_eV(), obj.get_weights());
    }

    static boost::python::tuple getstate(boost::python::object obj) {
      return boost::python::make_tuple(obj.attr("__dict__"));
    }

    static void setstate(boost::python::object obj, boost::python::tuple state) {
      DXTBX_ASSERT(boost::python::len(state) == 2);

      // restore the object's __dict__
      boost::python::dict d =
        boost::python::extract<boost::python::dict>(obj.attr("__dict__"))();
      d.update(state[0]);
    }

    static bool getstate_manages_dict() {
      return true;
    }
  };

  template <>
  boost::python::dict to_dict<Spectrum>(const Spectrum &obj) {
    boost::python::dict result;
    result["energies"] = obj.get_energies_eV();
    result["weights"] = obj.get_weights();
    return result;
  }

  template <>
  Spectrum *from_dict<Spectrum>(boost::python::dict obj) {
    Spectrum *s = new Spectrum(boost::python::extract<vecd>(obj["energies"]),
                               boost::python::extract<vecd>(obj["weights"]));
    return s;
  }

  void export_spectrum() {
    // Export Spectrum
    class_<Spectrum, boost::shared_ptr<Spectrum> >("Spectrum")
      .def(init<const Spectrum &>())
      .def(init<vecd, vecd>((arg("energies"), arg("weights"))))
      .def("get_energies_eV", &Spectrum::get_energies_eV)
      .def("get_weights", &Spectrum::get_weights)
      .def("get_weighted_energy_eV", &Spectrum::get_weighted_energy_eV)
      .def("get_weighted_wavelength", &Spectrum::get_weighted_wavelength)
      .def("get_emin_ev", &Spectrum::get_emin_ev)
      .def("get_emax_ev", &Spectrum::get_emax_ev)
      .def("__str__", &spectrum_to_string)
      .def("to_dict", &to_dict<Spectrum>)
      .def("from_dict", &from_dict<Spectrum>, return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(SpectrumPickleSuite());

    scitbx::af::boost_python::flex_wrapper<Spectrum>::plain("flex_Spectrum");
  }

}}}  // namespace dxtbx::model::boost_python
