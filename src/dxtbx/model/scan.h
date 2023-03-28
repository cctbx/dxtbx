/*
 * scan.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_SCAN_H
#define DXTBX_MODEL_SCAN_H

#include <cmath>
#include <iostream>
#include <map>
#include <scitbx/vec2.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include "scan_helpers.h"
#include <dxtbx/array_family/flex_table.h>
#include <dxtbx/array_family/flex_table_suite.h>
#include <algorithm>

namespace dxtbx { namespace model {

  using dxtbx::af::flex_table;
  using dxtbx::af::flex_type_generator;
  using scitbx::mat3;
  using scitbx::rad_as_deg;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::constants::pi;

  typedef std::map<std::string, scitbx::af::shared<vec2<int> > > ExpImgRangeMap;

  typedef flex_type_generator<bool,
                              int,
                              std::size_t,
                              double,
                              std::string,
                              vec2<double>,
                              vec3<double>,
                              mat3<double>,
                              int6>::type scan_property_types;

  typedef dxtbx::af::flex_table<scan_property_types>::const_iterator const_iterator;

  class ScanBase {};

  class Scan : public ScanBase {
  public:
    Scan()
        : image_range_(0, 0),
          num_images_(1),
          batch_offset_(0),
          properties_(flex_table<scan_property_types>(1)) {}

    /**
     * @param image_range The range of images covered by the scan
     * @param batch_offset An offset to add to the image number (for tracking of
     *                     unique batch numbers for multi-crystal datasets)
     */
    Scan(vec2<int> image_range, int batch_offset = 0)
        : image_range_(image_range),
          num_images_(1 + image_range_[1] - image_range_[0]),
          batch_offset_(batch_offset) {
      DXTBX_ASSERT(num_images_ >= 0);
      properties_ = flex_table<scan_property_types>(num_images_);
    }

    /**
     * @param image_range The range of images covered by the scan
     * @param oscillation A tuple containing the start angle of the first image
     *                    and the oscillation range (the angular width) of each
     *                    frame
     * @param batch_offset An offset to add to the image number (for tracking of
     *                     unique batch numbers for multi-crystal datasets)
     */
    Scan(vec2<int> image_range, vec2<double> oscillation, int batch_offset = 0)
        : image_range_(image_range),
          num_images_(1 + image_range_[1] - image_range_[0]),
          batch_offset_(batch_offset) {
      DXTBX_ASSERT(num_images_ >= 0);
      properties_ = flex_table<scan_property_types>(num_images_);
      scitbx::af::shared<double> exposure_times =
        scitbx::af::shared<double>(num_images_, 0.0);
      set_property("exposure_time", exposure_times.const_ref());
      scitbx::af::shared<double> epochs = scitbx::af::shared<double>(num_images_, 0.0);
      set_property("epochs", epochs.const_ref());
      set_oscillation(oscillation);
      DXTBX_ASSERT(properties_.is_consistent());
    }

    /**
     * @param image_range The range of images covered by the scan
     * @param oscillation A tuple containing the start angle of the first image
     *                    and the oscillation range (the angular width) of each
     *                    frame
     * @param exposure_times The exposure duration of each image
     * @param epochs The time of capture for each image
     * @param batch_offset A offset to add to the image number (for tracking of
     *                     unique batch numbers for multi-crystal datasets)
     */
    Scan(vec2<int> image_range,
         vec2<double> oscillation,
         const scitbx::af::shared<double> &exposure_times,
         const scitbx::af::shared<double> &epochs,
         int batch_offset = 0)
        : image_range_(image_range),
          num_images_(1 + image_range_[1] - image_range_[0]),
          batch_offset_(batch_offset) {
      DXTBX_ASSERT(num_images_ >= 0);
      DXTBX_ASSERT(oscillation[1] >= 0);
      properties_ = flex_table<scan_property_types>(num_images_);

      if (exposure_times.size() == 1 && num_images_ > 1) {
        // assume same exposure time for all images - there is
        // probably a better way of coding this...
        scitbx::af::shared<double> expanded_exposure_times;
        expanded_exposure_times.reserve(num_images_);
        for (int j = 0; j < num_images_; j++) {
          expanded_exposure_times.push_back(exposure_times[0]);
        }
        set_property("exposure_time", expanded_exposure_times.const_ref());
      } else {
        set_property("exposure_time", exposure_times.const_ref());
      }
      set_property("epochs", epochs.const_ref());
      set_oscillation(oscillation);
      DXTBX_ASSERT(properties_.is_consistent());
    }

    /**
     * @param image_range The range of images covered by the scan
     * @param properties_table Hash table of different properties for each image
     * @param batch_offset A offset to add to the image number (for tracking of
     *                     unique batch numbers for multi-crystal datasets)
     */
    Scan(vec2<int> image_range,
         flex_table<scan_property_types> properties_table,
         int batch_offset = 0)
        : image_range_(image_range),
          num_images_(1 + image_range_[1] - image_range_[0]),
          batch_offset_(batch_offset) {
      DXTBX_ASSERT(num_images_ >= 0);
      DXTBX_ASSERT(properties_table.is_consistent());
      DXTBX_ASSERT(properties_table.size() == num_images_);
      properties_ = properties_table;
    }

    /** Copy */
    Scan(const Scan &rhs)
        : image_range_(rhs.image_range_),
          valid_image_ranges_(rhs.valid_image_ranges_),
          num_images_(rhs.num_images_),
          batch_offset_(rhs.batch_offset_),
          properties_(rhs.properties_) {}

    virtual ~Scan() {}

    vec2<int> get_image_range() const {
      return image_range_;
    }

    /** Not exported to python **/
    ExpImgRangeMap get_valid_image_ranges_map() const {
      return valid_image_ranges_;
    }

    /** Get the element for a given key if it exists, else return empty array**/
    scitbx::af::shared<vec2<int> > get_valid_image_ranges_key(std::string i) const {
      typedef ExpImgRangeMap::const_iterator iterator;
      for (iterator it = valid_image_ranges_.begin(); it != valid_image_ranges_.end();
           ++it) {
        if (it->first == i) {
          return it->second;
        }
      }
      scitbx::af::shared<vec2<int> > empty;
      return empty;
    }

    void set_valid_image_ranges_array(std::string i,
                                      scitbx::af::shared<vec2<int> > values) {
      /** Set a list of valid image range tuples for experiment identifier 'i'**/
      for (std::size_t j = 0; j < values.size(); ++j) {
        vec2<int> pair = values[j];
        DXTBX_ASSERT(pair[0] >= image_range_[0]);
        DXTBX_ASSERT(pair[0] <= image_range_[1]);
        DXTBX_ASSERT(pair[1] >= image_range_[0]);
        DXTBX_ASSERT(pair[1] <= image_range_[1]);
      }
      valid_image_ranges_[i] = values;
    }

    int get_batch_offset() const {
      return batch_offset_;
    }

    bool contains(std::string property_name) const {
      return properties_.contains(property_name);
    }

    template <typename T>
    scitbx::af::shared<T> get_property(
      const flex_table<scan_property_types>::key_type &key) const {
      DXTBX_ASSERT(properties_.contains(key));
      return properties_.get<T>(key);
    }

    template <typename T>
    void set_property(const typename flex_table<scan_property_types>::key_type &key,
                      const scitbx::af::const_ref<T> &value) {
      DXTBX_ASSERT(value.size() == properties_.size());
      dxtbx::af::flex_table_suite::setitem_column(properties_, key, value);
    }

    flex_table<scan_property_types> get_properties() const {
      return properties_;
    }

    void set_properties(flex_table<scan_property_types> new_table) {
      DXTBX_ASSERT(new_table.is_consistent());
      DXTBX_ASSERT(new_table.size() == num_images_);
      properties_ = new_table;
    }

    bool is_still() const {
      if (!properties_.contains("oscillation")) {
        return false;
      }
      if (properties_.size() == 0) {
        return false;
      }

      return std::abs(properties_.get<double>("oscillation")[1])
             < min_oscillation_width_;
    }

    int get_batch_for_image_index(int index) const {
      return index + batch_offset_;
    }

    int get_batch_for_array_index(int index) const {
      return index + batch_offset_ + 1;
    }

    vec2<int> get_batch_range() const {
      return vec2<int>(image_range_[0] + batch_offset_,
                       image_range_[1] + batch_offset_);
    }

    /** (zero based) */
    vec2<int> get_array_range() const {
      return vec2<int>(image_range_[0] - 1, image_range_[1]);
    }

    vec2<double> get_oscillation() const {
      DXTBX_ASSERT(properties_.contains("oscillation"));
      scitbx::af::shared<double> osc = properties_.get<double>("oscillation");

      // Special case for Scans with single image
      if (properties_.size() == 1) {
        DXTBX_ASSERT(properties_.contains("oscillation_width"));
        scitbx::af::shared<double> osc_width =
          properties_.get<double>("oscillation_width");
        return vec2<double>(osc[0], osc_width[0]);
      }

      DXTBX_ASSERT(properties_.size() > 1);
      return vec2<double>(osc[0], osc[1]);
    }

    vec2<double> get_oscillation_in_deg() const {
      DXTBX_ASSERT(properties_.contains("oscillation"));
      scitbx::af::shared<double> osc = properties_.get<double>("oscillation");
      return vec2<double>(rad_as_deg(osc[0]), rad_as_deg(osc[1]));
    }

    scitbx::af::shared<double> get_oscillation_arr_in_deg() const {
      DXTBX_ASSERT(properties_.contains("oscillation"));
      scitbx::af::shared<double> osc = properties_.get<double>("oscillation");
      scitbx::af::shared<double> osc_in_deg;
      osc_in_deg.resize(osc.size());
      std::transform(osc.begin(), osc.end(), osc_in_deg.begin(), rad_as_deg);
      return osc_in_deg;
    }

    int get_num_images() const {
      return num_images_;
    }

    scitbx::af::shared<double> get_exposure_times() const {
      DXTBX_ASSERT(properties_.contains("exposure_time"));
      return properties_.get<double>("exposure_time");
    }

    scitbx::af::shared<double> get_epochs() const {
      DXTBX_ASSERT(properties_.contains("epochs"));
      return properties_.get<double>("epochs");
    }

    void set_image_range(vec2<int> image_range) {
      /*
       * Known issue with this function
       * https://github.com/cctbx/dxtbx/issues/497
       */
      image_range_ = image_range;
      num_images_ = 1 + image_range_[1] - image_range_[0];
      properties_.resize(num_images_);
      DXTBX_ASSERT(num_images_ > 0);
    }

    void set_batch_offset(int batch_offset) {
      batch_offset_ = batch_offset;
    }

    void set_oscillation(vec2<double> oscillation) {
      DXTBX_ASSERT(oscillation[1] >= 0.0);
      scitbx::af::shared<double> oscillation_arr(num_images_);
      for (std::size_t i = 0; i < num_images_; ++i) {
        oscillation_arr[i] = oscillation[0] + oscillation[1] * i;
      }
      set_property("oscillation", oscillation_arr.const_ref());

      // Special case for Scan for one image
      if (num_images_ == 1) {
        scitbx::af::shared<double> oscillation_width_arr(num_images_);
        oscillation_width_arr[0] = oscillation[1];
        set_property("oscillation_width", oscillation_width_arr.const_ref());
      }
    }

    void set_exposure_times(scitbx::af::shared<double> exposure_times) {
      DXTBX_ASSERT(exposure_times.size() == num_images_);
      set_property("exposure_time", exposure_times.const_ref());
      DXTBX_ASSERT(properties_.is_consistent());
    }

    void set_epochs(const scitbx::af::shared<double> &epochs) {
      DXTBX_ASSERT(epochs.size() == num_images_);
      set_property("epochs", epochs.const_ref());
      DXTBX_ASSERT(properties_.is_consistent());
    }

    /** Get the total oscillation range of the scan */
    vec2<double> get_oscillation_range() const {
      vec2<double> osc = get_oscillation();
      return vec2<double>(osc[0], osc[0] + num_images_ * osc[1]);
    }

    /** Get the image angle and oscillation width as a tuple */
    vec2<double> get_image_oscillation(int index) const {
      DXTBX_ASSERT(image_range_[0] <= index && index <= image_range_[1]);
      vec2<double> osc = get_oscillation();
      return vec2<double>(osc[0] + (index - image_range_[0]) * osc[1], osc[1]);
    }

    double get_image_epoch(int index) const {
      DXTBX_ASSERT(properties_.contains("epochs"));
      DXTBX_ASSERT(image_range_[0] <= index && index <= image_range_[1]);
      return properties_.get<double>("epochs")[index - image_range_[0]];
    }

    double get_image_exposure_time(int index) const {
      DXTBX_ASSERT(properties_.contains("exposure_time"));
      DXTBX_ASSERT(image_range_[0] <= index && index <= image_range_[1]);
      return properties_.get<double>("exposure_time")[index - image_range_[0]];
    }

    bool operator==(const Scan &rhs) const {
      if (image_range_ != image_range_ || batch_offset_ != rhs.batch_offset_) {
        return false;
      }

      if (properties_.size() != rhs.properties_.size()) {
        return false;
      }
      for (const_iterator it = properties_.begin(); it != properties_.end(); ++it) {
        if (!rhs.properties_.contains(it->first)) {
          return false;
        }
        // TODO explicit column comparison
      }
      return true;
    }

    bool operator!=(const Scan &scan) const {
      return !(*this == scan);
    }

    bool operator<(const Scan &scan) const {
      return image_range_[0] < scan.image_range_[0];
    }

    bool operator<=(const Scan &scan) const {
      return image_range_[0] <= scan.image_range_[0];
    }

    bool operator>(const Scan &scan) const {
      return image_range_[0] > scan.image_range_[0];
    }

    bool operator>=(const Scan &scan) const {
      return image_range_[0] >= scan.image_range_[0];
    }

    /**
     * Append the rhs scan onto the current scan
     */

    void append(const Scan &rhs, double scan_tolerance) {
      using namespace dxtbx::af::flex_table_suite;
      DXTBX_ASSERT(image_range_[1] + 1 == rhs.image_range_[0]);
      DXTBX_ASSERT(batch_offset_ == rhs.batch_offset_);

      flex_table<scan_property_types> rhs_properties = rhs.get_properties();

      // Explicitly check oscillation
      if (properties_.contains("oscillation")) {
        double osc_width = get_oscillation()[1];
        double eps = scan_tolerance * std::abs(osc_width);
        double rhs_osc_width = rhs.get_oscillation()[1];
        DXTBX_ASSERT(eps > 0);
        DXTBX_ASSERT(std::abs(osc_width) > min_oscillation_width_);
        DXTBX_ASSERT(std::abs(osc_width - rhs_osc_width) < eps);
        // sometimes ticking through 0 the first difference is not helpful
        double diff_2pi = std::abs(mod_2pi(get_oscillation_range()[1])
                                   - mod_2pi(rhs.get_oscillation_range()[0]));
        double diff_abs =
          std::abs(get_oscillation_range()[1] - rhs.get_oscillation_range()[0]);
        DXTBX_ASSERT(std::min(diff_2pi, diff_abs) < eps * get_num_images());

        /*
        If either properties table contains oscillation_width, remove as this
        is only needed for scans of one image
        */
        if (properties_.contains("oscillation_width")) {
          delitem_column(properties_, "oscillation_width");
        }
        if (rhs_properties.contains("oscillation_width")) {
          delitem_column(rhs_properties, "oscillation_width");
        }
      }

      image_range_[1] = rhs.image_range_[1];
      num_images_ = 1 + image_range_[1] - image_range_[0];

      for (const_iterator it = properties_.begin(); it != properties_.end(); ++it) {
        DXTBX_ASSERT(rhs_properties.contains(it->first));
      }
      extend(properties_, rhs_properties);
    }

    /**
     * Append the rhs scan onto the current scan
     */
    Scan &operator+=(const Scan &rhs) {
      // Set the epsilon to 1% of oscillation range
      append(rhs, 0.01);
      return *this;
    }

    /**
     * Return a new scan which consists of the contents of this scan and
     * the contents of the other scan, provided that they are consistent.
     * If they are not consistent then an AssertionError will result.
     */
    Scan operator+(const Scan &rhs) const {
      Scan lhs(*this);
      lhs += rhs;
      return lhs;
    }

    /**
     * Check if the angle is in the range of angles covered by the scan.
     */
    bool is_angle_valid(double angle) const {
      return is_angle_in_range(get_oscillation_range(), angle);
    }

    bool is_image_index_valid(double index) const {
      return (image_range_[0] <= index && index <= image_range_[1]);
    }

    bool is_batch_valid(int batch) const {
      vec2<int> batch_range = get_batch_range();
      return (batch_range[0] <= batch && batch <= batch_range[1]);
    }

    bool is_array_index_valid(double index) const {
      return is_image_index_valid(index + 1);
    }

    /**
     * Calculate the angle corresponding to the given frame
     * @param index The frame number
     * @returns The angle at the given frame
     */
    double get_angle_from_image_index(double index) const {
      vec2<double> oscillation = get_oscillation();
      return oscillation[0] + (index - image_range_[0]) * oscillation[1];
    }

    /**
     * Calculate the angle corresponding to the given zero based frame
     * @param index The frame number
     * @returns The angle at the given frame
     */
    double get_angle_from_array_index(double index) const {
      return get_angle_from_image_index(index + 1);
    }

    /**
     * Calculate the frame corresponding to the given angle
     * @param angle The angle
     * @returns The frame at the given angle
     */
    double get_image_index_from_angle(double angle) const {
      vec2<double> oscillation = get_oscillation();
      return image_range_[0] + (angle - oscillation[0]) / oscillation[1];
    }

    /**
     * Calculate the zero based frame corresponding to the given angle
     * @param angle The angle
     * @returns The frame at the given angle
     */
    double get_array_index_from_angle(double angle) const {
      return get_image_index_from_angle(angle) - 1;
    }

    /**
     * A function to calculate all the frames in the scan at which an
     * observation with a given angle will be observed. I.e. for a given angle,
     * find all the equivalent angles (i.e. mod 2pi) within the scan range and
     * calculate the frame number for each angle.
     * Calculate and return an array of frame numbers at which a reflection
     * with a given rotation angle will be observed.
     * @param angle The rotation angle of the reflection
     * @returns The array of frame numbers
     */
    scitbx::af::shared<vec2<double> > get_image_indices_with_angle(double angle) const {
      scitbx::af::shared<double> angles =
        get_mod2pi_angles_in_range(get_oscillation_range(), angle);
      scitbx::af::shared<vec2<double> > result(angles.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i][0] = angles[i];
        result[i][1] = get_image_index_from_angle(angles[i]);
      }
      return result;
    }

    /**
     * Calculate and return an array of zero based frame numbers at which a
     * reflection with a given rotation angle will be observed.
     * @param angle The rotation angle of the reflection
     * @returns The array of frame numbers
     */
    scitbx::af::shared<vec2<double> > get_array_indices_with_angle(
      double angle,
      double padding = 0,
      bool deg = false) const {
      DXTBX_ASSERT(padding >= 0);
      if (deg == true) {
        padding = padding * pi / 180.0;
      }
      vec2<double> range = get_oscillation_range();
      range[0] -= padding;
      range[1] += padding;
      scitbx::af::shared<double> angles = get_mod2pi_angles_in_range(range, angle);
      scitbx::af::shared<vec2<double> > result(angles.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i][0] = angles[i];
        result[i][1] = get_array_index_from_angle(angles[i]);
      }
      return result;
    }

    Scan operator[](int index) const {
      using namespace dxtbx::af::flex_table_suite;
      DXTBX_ASSERT((index >= 0) && (index < get_num_images()));

      int image_index = image_range_[0] + index;
      flex_table<scan_property_types> properties_slice =
        get_properties_slice(boost::python::slice(index, index + 1));

      return Scan(
        vec2<int>(image_index, image_index), properties_slice, get_batch_offset());
    }

    flex_table<scan_property_types> get_properties_slice(
      boost::python::slice slice) const {
      /*
        Wrapper for slicing the properties table to account for oscillation
        width being lost in the slice if a Scan with a single image
      */

      if (properties_.contains("oscillation")) {
        double oscillation_width = get_oscillation()[1];
        flex_table<scan_property_types> sliced_properties =
          dxtbx::af::flex_table_suite::getitem_slice(properties_, slice);
        if (sliced_properties.size() == 1) {
          scitbx::af::shared<double> oscillation_width_arr(1);
          oscillation_width_arr[0] = oscillation_width;
          dxtbx::af::flex_table_suite::setitem_column(
            sliced_properties, "oscillation_width", oscillation_width_arr.const_ref());
          return sliced_properties;
        }
        return sliced_properties;
      }
      return dxtbx::af::flex_table_suite::getitem_slice(properties_, slice);
    }

    friend std::ostream &operator<<(std::ostream &os, const Scan &s);

  private:
    vec2<int> image_range_;
    ExpImgRangeMap valid_image_ranges_; /** initialised as an empty map **/
    double min_oscillation_width_ = 1e-7;
    int num_images_;
    int batch_offset_;
    flex_table<scan_property_types> properties_;
  };

  /*
   * Summary for operator<<
  void add_property_summary(std::ostream &os,
                            const std::string property_name,
                            const flex_table<scan_property_types> &properties,
                            Scan &s) {
    if (property_name == "oscillation") {
      DXTBX_ASSERT(properties.contains("oscillation"));
      vec2<double> oscillation = s.get_oscillation_in_deg();
      os << "    oscillation:   " << oscillation.const_ref() << "\n";
      return;
    } else if (property_name == "exposure_time") {
      DXTBX_ASSERT(properties.contains("exposure_time"));
      os << "    exposure time: "
         << properties.get<double>("exposure_time").const_ref()[0] << "\n";
      return;
    } else if (property_name == "epochs") {
      DXTBX_ASSERT(properties.contains("epochs"));
      os << "    init epoch: " << properties.get<double>("epochs").const_ref()[0]
         << "\n";
      return;
    }
    DXTBX_ERROR("No summary found for " + (property_name));
  }
  */
  /** Print Scan information */
  inline std::ostream &operator<<(std::ostream &os, const Scan &s) {
    os << "Scan:\n";
    os << "    number of images:   " << s.get_num_images() << "\n";
    os << "    image range:   " << s.get_image_range().const_ref() << "\n";
    /*
    flex_table<scan_property_types> properties = s.get_properties();
    for (const_iterator it = properties.begin(); it != properties.end(); ++it) {
      add_property_summary(os, it->first, properties, s);
    }
    */
    return os;
  }

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_SCAN_H
