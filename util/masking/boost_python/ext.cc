#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dxtbx/util/masking/masking.h>
#include <dxtbx/util/masking/goniometer_shadow_masking.h>


namespace dxtbx { namespace util { namespace masking { namespace boost_python {

  using scitbx::vec2;
  using scitbx::vec3;

  scitbx::af::shared<bool> is_inside_polygon_a(const scitbx::af::const_ref<vec2<double> > &poly,
                                       const scitbx::af::const_ref<vec2<double> > &xy) {
    scitbx::af::shared<bool> inside(xy.size(), false);
    for (std::size_t i = 0; i < xy.size(); i++) {
      inside[i] = is_inside_polygon(poly, xy[i][0], xy[i][1]);
    }
    return inside;
  }

  static boost::python::list GoniometerShadowMasker_project_extrema(
    GoniometerShadowMasker &masker,
    const Detector &detector,
    double scan_angle) {
    scitbx::af::shared<scitbx::af::shared<scitbx::vec2<double> > > result =
      masker.project_extrema(detector, scan_angle);
    typename scitbx::af::shared<scitbx::af::shared<scitbx::vec2<double> > >::iterator
      iter;
    boost::python::list list;
    for (iter = result.begin(); iter != result.end(); ++iter) {
      list.append(*iter);
    }
    return list;
  }

  // Copied from dxtbx/boost_python/imageset_ext.cc
  namespace detail {
    /**
     * Tuple from list
     */
    boost::python::tuple list_to_tuple(boost::python::list x) {
      boost::python::object main = boost::python::import("__main__");
      boost::python::object global(main.attr("__dict__"));
      boost::python::object result =
        exec("def func_list_to_tuple(x):return tuple(x)", global, global);
      boost::python::object func_list_to_tuple = global["func_list_to_tuple"];
      return boost::python::extract<boost::python::tuple>(func_list_to_tuple(x))();
    }
  }  // namespace detail

  // Copied from dxtbx/boost_python/imageset_ext.cc
  template <typename T>
  boost::python::tuple image_as_tuple(const Image<T> &image) {
    boost::python::list result;
    for (std::size_t i = 0; i < image.n_tiles(); ++i) {
      result.append(image.tile(i).data());
    }
    return detail::list_to_tuple(result);
  }

  boost::python::tuple GoniometerShadowMasker_get_mask(

    GoniometerShadowMasker &masker,
    const Detector &detector,
    double scan_angle) {
    return image_as_tuple<bool>(masker.get_mask(detector, scan_angle));
  }

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dxtbx_util_masking_ext) {

    def("mask_untrusted_rectangle", &mask_untrusted_rectangle);

    def("mask_untrusted_circle", &mask_untrusted_circle);

    def("mask_untrusted_resolution_range", &mask_untrusted_resolution_range);

    def("mask_untrusted_polygon", &mask_untrusted_polygon);

    def("is_inside_polygon", &is_inside_polygon);

    def("is_inside_polygon", &is_inside_polygon_a);

    class_<GoniometerShadowMasker>("GoniometerShadowMasker", no_init)
      .def(init<const MultiAxisGoniometer &,
                const scitbx::af::const_ref<scitbx::vec3<double> > &,
                const scitbx::af::const_ref<std::size_t> &>())
      .def("extrema_at_scan_angle",
           &GoniometerShadowMasker::extrema_at_scan_angle)
      .def("set_goniometer_angles",
           &GoniometerShadowMasker::set_goniometer_angles)
      .def("project_extrema", GoniometerShadowMasker_project_extrema)
      .def("get_mask", GoniometerShadowMasker_get_mask);

    class_<SmarGonShadowMasker, bases<GoniometerShadowMasker> >(
      "SmarGonShadowMasker", no_init)
      .def(init<const MultiAxisGoniometer &>())
      .def("extrema_at_scan_angle", &SmarGonShadowMasker::extrema_at_scan_angle);

  }
}}}}  // namespace dxtbx::util::masking::boost_python
