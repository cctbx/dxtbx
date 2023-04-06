/*
 * scan.cc
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
#include <boost/python/make_constructor.hpp>
#include <boost/python/slice.hpp>
#include <memory>
#include <string>
#include <sstream>
#include <scitbx/constants.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/boost_python/to_from_dict.h>
#include <dxtbx/array_family/flex_table.h>
#include <dxtbx/array_family/flex_table_suite.h>
#include <algorithm>

namespace dxtbx { namespace model { namespace boost_python {

  using dxtbx::model::scan_property_types;
  using namespace boost::python;
  using namespace dxtbx::af::flex_table_suite;
  using scitbx::deg_as_rad;
  using scitbx::rad_as_deg;

  static vec2<double> rad_as_deg(vec2<double> angles) {
    angles[0] = rad_as_deg(angles[0]);
    angles[1] = rad_as_deg(angles[1]);
    return angles;
  }

  static vec2<double> deg_as_rad(vec2<double> angles) {
    angles[0] = deg_as_rad(angles[0]);
    angles[1] = deg_as_rad(angles[1]);
    return angles;
  }

  static scitbx::af::shared<double> deg_as_rad(
    scitbx::af::shared<double> angles_in_deg) {
    scitbx::af::shared<double> angles_in_rad;
    angles_in_rad.resize(angles_in_deg.size());
    std::transform(angles_in_deg.begin(),
                   angles_in_deg.end(),
                   angles_in_rad.begin(),
                   scitbx::deg_as_rad);
    return angles_in_rad;
  }

  std::string scan_to_string(const Scan &scan) {
    std::stringstream ss;
    ss << scan;
    return ss.str();
  }

  flex_table<scan_property_types> extract_properties_table(
    boost::python::dict properties_dict,
    int num_images,
    bool convert_oscillation_to_rad = false) {
    boost::python::list keys = properties_dict.keys();
    boost::python::list values = boost::python::list(properties_dict.values());
    flex_table<scan_property_types> properties =
      flex_table<scan_property_types>(num_images);

    for (int i = 0; i < len(keys); ++i) {
      std::string key = boost::python::extract<std::string>(keys[i]);
      object value = values[i];
      DXTBX_ASSERT(len(value) == num_images);
      std::string obj_type = boost::python::extract<std::string>(
        value[0].attr("__class__").attr("__name__"));

      // Handled explicitly as it is in deg when serialised but rad in code
      if (key == "oscillation" && convert_oscillation_to_rad) {
        DXTBX_ASSERT(obj_type == "float");
        scitbx::af::shared<double> osc_in_deg =
          boost::python::extract<scitbx::af::shared<double> >(value);
        properties[key] = deg_as_rad(osc_in_deg);
      } else if (obj_type == "int") {
        properties[key] = boost::python::extract<scitbx::af::shared<int> >(value);
      } else if (obj_type == "float") {
        properties[key] = boost::python::extract<scitbx::af::shared<double> >(value);
      } else if (obj_type == "bool") {
        properties[key] = boost::python::extract<scitbx::af::shared<bool> >(value);
      } else if (obj_type == "str") {
        properties[key] =
          boost::python::extract<scitbx::af::shared<std::string> >(value);
      } else if (obj_type == "tuple") {
        std::string element_type = boost::python::extract<std::string>(
          value[0][0].attr("__class__").attr("__name__"));
        int tuple_size = len(value[0]);

        if (tuple_size == 2) {
          if (element_type == "float") {
            properties[key] =
              boost::python::extract<scitbx::af::shared<vec2<double> > >(value);
          } else {
            throw DXTBX_ERROR("Unknown type for column name " + key);
          }
        } else if (tuple_size == 3) {
          if (element_type == "float") {
            properties[key] =
              boost::python::extract<scitbx::af::shared<vec3<double> > >(value);
          } else {
            throw DXTBX_ERROR("Unknown type for column name " + key);
          }
        } else {
          throw DXTBX_ERROR("Unknown type for column name " + key);
        }
      } else {
        throw DXTBX_ERROR("Unknown type for column name " + key);
      }
    }

    return properties;
  }

  void set_properties_table_from_dict(Scan &obj, boost::python::dict properties_dict) {
    int num_images = len(boost::python::list(properties_dict.values())[0]);
    flex_table<scan_property_types> properties =
      extract_properties_table(properties_dict, num_images);
    obj.set_properties(properties);
  }

  struct ScanPickleSuite : boost::python::pickle_suite {
    typedef flex_table<scan_property_types>::const_iterator const_iterator;

    static boost::python::tuple getinitargs(const Scan &obj) {
      return boost::python::make_tuple(obj.get_image_range(), obj.get_batch_offset());
    }
    static boost::python::tuple getstate(const Scan &obj) {
      flex_table<scan_property_types> properties = obj.get_properties();
      boost::python::dict properties_dict;
      column_to_object_visitor visitor;
      for (const_iterator it = properties.begin(); it != properties.end(); ++it) {
        properties_dict[it->first] = it->second.apply_visitor(visitor);
      }

      return boost::python::make_tuple(
        properties.nrows(), properties.ncols(), properties_dict);
    }

    static void setstate(Scan &obj, boost::python::tuple state) {
      DXTBX_ASSERT(boost::python::len(state) == 3);
      std::size_t nrows = extract<std::size_t>(state[0]);
      std::size_t ncols = extract<std::size_t>(state[1]);
      boost::python::dict properties_dict = extract<dict>(state[2]);
      DXTBX_ASSERT(len(properties_dict) == ncols);
      flex_table<scan_property_types> properties =
        extract_properties_table(properties_dict, nrows);
      obj.set_properties(properties);
    }
  };

  template <typename T>
  struct scan_property_table_wrapper
      : public dxtbx::af::flex_table_suite::flex_table_wrapper<T> {
    typedef dxtbx::af::flex_table_suite::flex_table_wrapper<T> base_type;
    typedef typename base_type::flex_table_type flex_table_type;
    typedef typename base_type::class_type class_type;

    static class_type wrap(const char *name) {
      return base_type::wrap(name);
    }
  };

  boost::python::dict MaptoPythonDict(ExpImgRangeMap map) {
    ExpImgRangeMap::iterator iter;
    boost::python::dict dictionary;
    for (iter = map.begin(); iter != map.end(); ++iter) {
      scitbx::af::shared<vec2<int> > val = iter->second;
      boost::python::list result;
      for (int k = 0; k < val.size(); ++k) {
        result.append(val[k]);
      }
      dictionary[iter->first] = result;
    }
    return dictionary;
  }

  boost::python::dict get_properties_dict(const Scan &obj) {
    boost::python::dict properties_dict;

    flex_table<scan_property_types> properties = obj.get_properties();
    column_to_object_visitor visitor;
    for (const_iterator it = properties.begin(); it != properties.end(); ++it) {
      properties_dict[it->first] =
        boost::python::tuple(it->second.apply_visitor(visitor));
    }
    return properties_dict;
  }

  template <>
  boost::python::dict to_dict<Scan>(const Scan &obj) {
    boost::python::dict result;
    result["image_range"] = obj.get_image_range();
    result["batch_offset"] = obj.get_batch_offset();

    flex_table<scan_property_types> properties = obj.get_properties();
    boost::python::dict properties_dict;
    column_to_object_visitor visitor;
    for (const_iterator it = properties.begin(); it != properties.end(); ++it) {
      if (it->first == "oscillation") {
        properties_dict[it->first] =
          boost::python::tuple(obj.get_oscillation_arr_in_deg());
      } else {
        properties_dict[it->first] =
          boost::python::tuple(it->second.apply_visitor(visitor));
      }
    }

    result["properties"] = properties_dict;

    boost::python::dict valid_image_ranges =
      MaptoPythonDict(obj.get_valid_image_ranges_map());
    result["valid_image_ranges"] = valid_image_ranges;
    return result;
  }

  inline scitbx::af::shared<double> make_exposure_times(std::size_t num,
                                                        boost::python::list obj) {
    scitbx::af::shared<double> result((scitbx::af::reserve(num)));
    std::size_t nl = boost::python::len(obj);
    DXTBX_ASSERT(num > 0 && nl <= num);
    if (nl == 0) {
      result.push_back(0.0);
      nl = 1;
    } else {
      for (std::size_t i = 0; i < nl; ++i) {
        result.push_back(boost::python::extract<double>(obj[i]));
      }
    }
    for (std::size_t i = nl; i < num; ++i) {
      result.push_back(result.back());
    }
    return result;
  }

  inline scitbx::af::shared<double> make_epochs(std::size_t num,
                                                boost::python::list obj) {
    scitbx::af::shared<double> result((scitbx::af::reserve(num)));
    std::size_t nl = boost::python::len(obj);
    DXTBX_ASSERT(num > 0 && nl <= num);
    if (nl == 0) {
      for (std::size_t i = 0; i < num; ++i) {
        result.push_back(0.0);
      }
    } else if (nl == 1) {
      for (std::size_t i = 0; i < num; ++i) {
        result.push_back(boost::python::extract<double>(obj[0]));
      }
    } else if (nl < num) {
      for (std::size_t i = 0; i < nl; ++i) {
        result.push_back(boost::python::extract<double>(obj[i]));
      }
      double e0 = result[result.size() - 1];
      double de = e0 - result[result.size() - 2];
      for (std::size_t i = 0; i < num - nl; ++i) {
        result.push_back(e0 + (i + 1) * de);
      }
    } else {
      for (std::size_t i = 0; i < num; ++i) {
        result.push_back(boost::python::extract<double>(obj[i]));
      }
    }
    return result;
  }

  scitbx::af::shared<double> make_oscillation_arr(std::size_t num_images,
                                                  vec2<double> oscillation) {
    scitbx::af::shared<double> oscillation_arr((scitbx::af::reserve(num_images)));
    for (std::size_t i = 0; i < num_images; ++i) {
      oscillation_arr.push_back(oscillation[0] + oscillation[1] * i);
    }
    return oscillation_arr;
  }

  template <>
  Scan *from_dict<Scan>(boost::python::dict obj) {
    vec2<int> ir = boost::python::extract<vec2<int> >(obj["image_range"]);
    int bo = boost::python::extract<int>(obj["batch_offset"]);
    DXTBX_ASSERT(ir[1] >= ir[0]);
    std::size_t num_images = ir[1] - ir[0] + 1;

    if (!obj.has_key("properties")) {
      /*
       * Legacy case with no properties table
       */

      vec2<double> osc =
        deg_as_rad(boost::python::extract<vec2<double> >(obj["oscillation"]));
      Scan *scan = new Scan(
        ir,
        osc,
        make_exposure_times(num_images,
                            boost::python::extract<boost::python::list>(
                              obj.get("exposure_time", boost::python::list()))),
        make_epochs(num_images,
                    boost::python::extract<boost::python::list>(
                      obj.get("epochs", boost::python::list()))),
        bo);
      boost::python::dict rangemap =
        boost::python::extract<boost::python::dict>(obj["valid_image_ranges"]);
      boost::python::list keys = rangemap.keys();
      boost::python::list values = rangemap.values();
      for (int i = 0; i < len(keys); ++i) {
        std::string key = boost::python::extract<std::string>(keys[i]);
        scitbx::af::shared<vec2<int> > result;
        int n_tuples = boost::python::len(values[i]);
        for (int n = 0; n < n_tuples; ++n) {
          result.push_back(boost::python::extract<vec2<int> >(values[i][n]));
        }
        scan->set_valid_image_ranges_array(key, result);
      }
      return scan;
    }

    Scan *scan = new Scan(ir, bo);

    boost::python::dict properties_dict =
      boost::python::extract<boost::python::dict>(obj["properties"]);
    scan->set_properties(extract_properties_table(properties_dict, num_images, true));

    boost::python::dict rangemap =
      boost::python::extract<boost::python::dict>(obj["valid_image_ranges"]);
    boost::python::list valid_img_keys = rangemap.keys();
    boost::python::list valid_img_values = rangemap.values();
    for (int i = 0; i < len(valid_img_keys); ++i) {
      std::string key = boost::python::extract<std::string>(valid_img_keys[i]);
      scitbx::af::shared<vec2<int> > result;
      int n_tuples = boost::python::len(valid_img_values[i]);
      for (int n = 0; n < n_tuples; ++n) {
        result.push_back(boost::python::extract<vec2<int> >(valid_img_values[i][n]));
      }
      scan->set_valid_image_ranges_array(key, result);
    }
    return scan;
  }

  static Scan scan_deepcopy(const Scan &scan, boost::python::object dict) {
    return Scan(scan);
  }

  static Scan scan_copy(const Scan &scan) {
    return Scan(scan);
  }

  static void set_valid_image_ranges(Scan &scan,
                                     std::string i,
                                     boost::python::list obj) {
    int n = boost::python::len(obj);
    scitbx::af::shared<vec2<int> > ranges;
    for (int k = 0; k < n; ++k) {
      ranges.push_back(boost_python::extract<vec2<int> >(obj[k]));
    }
    scan.set_valid_image_ranges_array(i, ranges);
  }

  static boost::python::list get_valid_image_ranges(Scan &scan, std::string i) {
    scitbx::af::shared<vec2<int> > ranges = scan.get_valid_image_ranges_key(i);
    boost::python::list result;
    if (ranges.size() != 0) {
      for (int k = 0; k < ranges.size(); ++k) {
        result.append(ranges[k]);
      }
    }
    return result;
  }

  static Scan *make_scan(vec2<int> image_range,
                         vec2<double> oscillation,
                         int batch_offset,
                         bool deg) {
    Scan *scan = NULL;
    if (deg) {
      scan =
        new Scan(image_range,
                 vec2<double>(deg_as_rad(oscillation[0]), deg_as_rad(oscillation[1])),
                 batch_offset);
    } else {
      scan = new Scan(image_range, oscillation, batch_offset);
    }
    return scan;
  }

  static Scan *make_scan_w_epoch(vec2<int> image_range,
                                 vec2<double> oscillation,
                                 const scitbx::af::shared<double> &exposure_times,
                                 const scitbx::af::shared<double> &epochs,
                                 int batch_offset,
                                 bool deg) {
    Scan *scan = NULL;
    if (deg) {
      scan =
        new Scan(image_range,
                 vec2<double>(deg_as_rad(oscillation[0]), deg_as_rad(oscillation[1])),
                 exposure_times,
                 epochs,
                 batch_offset);
    } else {
      scan = new Scan(image_range, oscillation, exposure_times, epochs, batch_offset);
    }
    return scan;
  }

  static Scan *make_scan_w_properties(vec2<int> image_range,
                                      boost::python::dict properties_dict,
                                      int batch_offset,
                                      bool deg) {
    int num_images = 1 + image_range[1] - image_range[0];
    flex_table<scan_property_types> properties =
      extract_properties_table(properties_dict, num_images);
    if (deg && properties.contains("oscillation")) {
      scitbx::af::shared<double> osc_in_deg = properties.get<double>("oscillation");
      scitbx::af::shared<double> osc = deg_as_rad(osc_in_deg);
      dxtbx::af::flex_table_suite::setitem_column(
        properties, "oscillation", osc.const_ref());
    }
    return new Scan(image_range, properties, batch_offset);
  }

  static Scan *make_scan_wo_properties(vec2<int> image_range, int batch_offset) {
    return new Scan(image_range, batch_offset);
  }

  static vec2<double> get_oscillation_range(const Scan &scan, bool deg) {
    vec2<double> range = scan.get_oscillation_range();
    return deg ? rad_as_deg(range) : range;
  }

  static vec2<double> get_oscillation(const Scan &scan, bool deg) {
    vec2<double> oscillation = scan.get_oscillation();
    return deg ? rad_as_deg(oscillation) : oscillation;
  }

  static void set_oscillation(Scan &scan, vec2<double> oscillation, bool deg) {
    if (deg) {
      oscillation = deg_as_rad(oscillation);
    }
    scan.set_oscillation(oscillation);
  }

  static vec2<double> get_image_oscillation(const Scan &scan, int image, bool deg) {
    vec2<double> oscillation = scan.get_image_oscillation(image);
    return deg ? rad_as_deg(oscillation) : oscillation;
  }

  static bool is_angle_valid(const Scan &scan, double angle, bool deg) {
    return scan.is_angle_valid(deg ? deg_as_rad(angle) : angle);
  }

  static scitbx::af::shared<bool> is_angle_valid_array(
    const Scan &scan,
    scitbx::af::const_ref<double> angle,
    bool deg) {
    scitbx::af::shared<bool> result(angle.size());
    for (std::size_t i = 0; i < angle.size(); ++i) {
      result[i] = scan.is_angle_valid(deg ? deg_as_rad(angle[i]) : angle[i]);
    }
    return result;
  }

  static double get_angle_from_image_index(const Scan &scan, double index, bool deg) {
    double angle = scan.get_angle_from_image_index(index);
    return deg ? rad_as_deg(angle) : angle;
  }

  static double get_angle_from_array_index(const Scan &scan, double index, bool deg) {
    double angle = scan.get_angle_from_array_index(index);
    return deg ? rad_as_deg(angle) : angle;
  }

  static scitbx::af::shared<double> get_angle_from_array_index_multiple(
    const Scan &scan,
    scitbx::af::const_ref<double> const &index,
    bool deg) {
    scitbx::af::shared<double> result((scitbx::af::reserve(index.size())));
    for (std::size_t i = 0; i < index.size(); i++) {
      result.push_back(get_angle_from_array_index(scan, index[i], deg));
    }
    return result;
  }

  static double get_image_index_from_angle(const Scan &scan, double angle, bool deg) {
    return scan.get_image_index_from_angle(deg ? deg_as_rad(angle) : angle);
  }

  static double get_array_index_from_angle(const Scan &scan, double angle, bool deg) {
    return scan.get_array_index_from_angle(deg ? deg_as_rad(angle) : angle);
  }

  static scitbx::af::shared<double> get_array_index_from_angle_multiple(
    const Scan &scan,
    scitbx::af::const_ref<double> const &angle,
    bool deg) {
    scitbx::af::shared<double> result((scitbx::af::reserve(angle.size())));
    for (std::size_t i = 0; i < angle.size(); i++) {
      result.push_back(get_array_index_from_angle(scan, angle[i], deg));
    }
    return result;
  }

  static scitbx::af::shared<vec2<double> >
  get_image_indices_with_angle(const Scan &scan, double angle, bool deg) {
    return scan.get_image_indices_with_angle(deg ? deg_as_rad(angle) : angle);
  }

  static scitbx::af::shared<vec2<double> >
  get_array_indices_with_angle(const Scan &scan, double angle, bool deg) {
    return scan.get_array_indices_with_angle(deg ? deg_as_rad(angle) : angle);
  }

  static Scan getitem_single(const Scan &scan, int index) {
    return scan[index];
  }

  static Scan getitem_slice(const Scan &scan, const slice index) {
    // Ensure no step
    DXTBX_ASSERT(index.step() == object());

    // Get start index
    int start = 0, stop = 0;
    if (index.start() == object()) {
      start = 0;
    } else {
      start = extract<int>(index.start());
    }

    // Get stop index
    if (index.stop() == object()) {
      stop = scan.get_num_images();
    } else {
      stop = extract<int>(index.stop());
    }

    // Check ranges
    DXTBX_ASSERT(start >= 0);
    DXTBX_ASSERT(stop <= scan.get_num_images());
    DXTBX_ASSERT(start < stop);

    // Create the new epoch array
    int first_image_index = scan.get_image_range()[0] + start;
    int last_image_index = scan.get_image_range()[0] + stop - 1;
    scitbx::af::shared<double> new_epochs(stop - start);
    for (std::size_t i = 0; i < new_epochs.size(); ++i) {
      new_epochs[i] = scan.get_image_epoch(first_image_index + i);
    }

    // Create the new epoch array
    scitbx::af::shared<double> new_exposure_times(stop - start);
    for (std::size_t i = 0; i < new_exposure_times.size(); ++i) {
      new_exposure_times[i] = scan.get_image_exposure_time(first_image_index + i);
    }

    // Create the new scan object
    return Scan(vec2<int>(first_image_index, last_image_index),
                scan.get_image_oscillation(first_image_index),
                new_exposure_times,
                new_epochs,
                scan.get_batch_offset());
  }

  void scan_swap(Scan &lhs, Scan &rhs) {
    std::swap(lhs, rhs);
  }

  template <typename T>
  object get_scan_property(const Scan &scan, const typename T::key_type &key) {
    DXTBX_ASSERT(scan.contains(key));
    flex_table<scan_property_types> properties = scan.get_properties();
    return getitem_column<T>(properties, key);
  }

  template <typename T>
  void set_scan_property(Scan &scan,
                         const std::string &key,
                         const scitbx::af::shared<T> &value) {
    scan.set_property(key, value.const_ref());
  }

  void export_scan() {
    scan_property_table_wrapper<flex_table<scan_property_types> >::wrap(
      "scan_property_table");

    // Export ScanBase
    class_<ScanBase>("ScanBase");

    // Export Scan : ScanBase
    class_<Scan, std::shared_ptr<Scan>, bases<ScanBase> >("Scan")
      .def(init<const Scan &>())
      .def("__init__",
           make_constructor(&make_scan,
                            default_call_policies(),
                            (arg("image_range"),
                             arg("oscillation"),
                             arg("batch_offset") = 0,
                             arg("deg") = true)))
      .def("__init__",
           make_constructor(&make_scan_w_epoch,
                            default_call_policies(),
                            (arg("image_range"),
                             arg("oscillation"),
                             arg("exposure_times"),
                             arg("epochs"),
                             arg("batch_offset") = 0,
                             arg("deg") = true)))
      .def("__init__",
           make_constructor(&make_scan_w_properties,
                            default_call_policies(),
                            (arg("image_range"),
                             arg("properties"),
                             arg("batch_offset") = 0,
                             arg("deg") = true)))
      .def("__init__",
           make_constructor(&make_scan_wo_properties,
                            default_call_policies(),
                            (arg("image_range"), arg("batch_offset") = 0)))
      .def("get_image_range", &Scan::get_image_range)
      .def("get_valid_image_ranges", get_valid_image_ranges)
      .def("set_valid_image_ranges", set_valid_image_ranges)
      .def("set_image_range", &Scan::set_image_range)
      .def("get_batch_offset", &Scan::get_batch_offset)
      .def("set_batch_offset", &Scan::set_batch_offset)
      .def("get_batch_for_image_index", &Scan::get_batch_for_image_index)
      .def("get_batch_for_array_index", &Scan::get_batch_for_array_index)
      .def("get_batch_range", &Scan::get_batch_range)
      .def("get_array_range", &Scan::get_array_range)
      .def("get_oscillation", &get_oscillation, (arg("deg") = true))
      .def("set_oscillation", &set_oscillation, (arg("deg") = true))
      .def("is_still", &Scan::is_still)
      .def("get_exposure_times", &Scan::get_exposure_times)
      .def("set_exposure_times", &Scan::set_exposure_times)
      .def("get_epochs", &Scan::get_epochs)
      .def("set_epochs", &Scan::set_epochs)
      .def("get_num_images", &Scan::get_num_images)
      .def("get_image_oscillation",
           &get_image_oscillation,
           (arg("index"), arg("deg") = true))
      .def("get_image_epoch", &Scan::get_image_epoch, (arg("index")))
      .def("get_oscillation_range", &get_oscillation_range, (arg("deg") = true))
      .def("is_angle_valid", &is_angle_valid, (arg("angle"), arg("deg") = true))
      .def("__deepcopy__", &scan_deepcopy)
      .def("__copy__", &scan_copy)
      .def("is_angle_valid", &is_angle_valid_array, (arg("angle"), arg("deg") = true))
      .def("is_image_index_valid", &Scan::is_image_index_valid, (arg("index")))
      .def("is_array_index_valid", &Scan::is_array_index_valid, (arg("index")))
      .def("is_batch_valid", &Scan::is_batch_valid, (arg("batch")))
      .def("get_angle_from_image_index",
           &get_angle_from_image_index,
           (arg("index"), arg("deg") = true))
      .def("get_angle_from_array_index",
           &get_angle_from_array_index,
           (arg("index"), arg("deg") = true))
      .def("get_angle_from_array_index",
           &get_angle_from_array_index_multiple,
           (arg("index"), arg("deg") = true))
      .def("get_image_index_from_angle",
           &get_image_index_from_angle,
           (arg("angle"), arg("deg") = true))
      .def("get_array_index_from_angle",
           &get_array_index_from_angle,
           (arg("angle"), arg("deg") = true))
      .def("get_array_index_from_angle",
           &get_array_index_from_angle_multiple,
           (arg("angle"), arg("deg") = true))
      .def("get_image_indices_with_angle",
           &get_image_indices_with_angle,
           (arg("angle"), arg("deg") = true))
      .def("get_array_indices_with_angle",
           &get_array_indices_with_angle,
           (arg("angle"), arg("deg") = true))
      .def("__getitem__", &getitem_single)
      .def("__getitem__", &getitem_slice)
      .def("get_property",
           &get_scan_property<flex_table<scan_property_types> >,
           (arg("key")))
      .def("get_properties", &get_properties_dict)
      .def("set_properties", &set_properties_table_from_dict, (arg("properties_dict")))
      .def("set_property", &set_scan_property<double>, (arg("key"), arg("value")))
      .def("set_property", &set_scan_property<bool>, (arg("key"), arg("value")))
      .def("set_property", &set_scan_property<int>, (arg("key"), arg("value")))
      .def("set_property", &set_scan_property<std::string>, (arg("key"), arg("value")))
      .def(
        "set_property", &set_scan_property<vec2<double> >, (arg("key"), arg("value")))
      .def(
        "set_property", &set_scan_property<vec3<double> >, (arg("key"), arg("value")))
      .def(self == self)
      .def(self != self)
      .def(self < self)
      .def(self <= self)
      .def(self > self)
      .def(self >= self)
      .def(self += self)
      .def(self + self)
      .def("append", &Scan::append, (arg("rhs"), arg("scan_tolerance") = 0.03))
      .def("__len__", &Scan::get_num_images)
      .def("__str__", &scan_to_string)
      .def("swap", &scan_swap)
      .def("to_dict", &to_dict<Scan>)
      .def("from_dict", &from_dict<Scan>, return_value_policy<manage_new_object>())
      .staticmethod("from_dict")
      .def_pickle(ScanPickleSuite());
  }

}}}  // namespace dxtbx::model::boost_python
