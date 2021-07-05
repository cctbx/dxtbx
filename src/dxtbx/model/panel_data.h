/*
 * panel_data.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_PANEL_DATA_H
#define DXTBX_MODEL_PANEL_DATA_H

#include <string>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/model_helpers.h>
#include <dxtbx/model/virtual_panel.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  using scitbx::af::int2;
  using scitbx::af::int4;
  using scitbx::af::tiny;

  /**
   * A panel class.
   */
  class PanelData : public VirtualPanel {
  public:
    /** Construct the panel with the simple px->mm strategy */
    PanelData()
        : pixel_size_(0.0, 0.0),
          image_size_(0, 0),
          trusted_range_(0.0, 0.0),
          thickness_(0.0),
          mu_(0.0),
          raw_image_offset_(0, 0) {}

    /** Construct with data */
    PanelData(std::string type,
              std::string name,
              tiny<double, 3> fast_axis,
              tiny<double, 3> slow_axis,
              tiny<double, 3> origin,
              tiny<double, 2> pixel_size,
              tiny<std::size_t, 2> image_size,
              tiny<double, 2> trusted_range,
              double thickness,
              std::string material,
              double mu)
        : pixel_size_(pixel_size),
          image_size_(image_size),
          trusted_range_(trusted_range),
          thickness_(thickness),
          material_(material),
          mu_(mu),
          raw_image_offset_(0, 0) {
      set_type(type);
      set_name(name);
      set_local_frame(fast_axis, slow_axis, origin);
    }

    virtual ~PanelData() {}

    /** Get the pixel size */
    tiny<double, 2> get_pixel_size() const {
      return pixel_size_;
    }

    /** Set the pixel size */
    void set_pixel_size(tiny<double, 2> pixel_size) {
      pixel_size_ = pixel_size;
    }

    /** Get the image size */
    tiny<std::size_t, 2> get_image_size() const {
      return image_size_;
    }

    /** Set the image size */
    void set_image_size(tiny<std::size_t, 2> image_size) {
      image_size_ = image_size;
    }

    /** Get the trusted range */
    tiny<double, 2> get_trusted_range() const {
      return trusted_range_;
    }

    /** Set the trusted range */
    void set_trusted_range(tiny<double, 2> trusted_range) {
      trusted_range_ = trusted_range;
    }

    /** Get the thickness */
    double get_thickness() const {
      return thickness_;
    }

    /** Set the thickness */
    void set_thickness(double thickness) {
      thickness_ = thickness;
    }

    /** Get the material */
    std::string get_material() const {
      return material_;
    }

    /** Set the material */
    void set_material(const std::string &material) {
      material_ = material;
    }

    /** Set mu */
    void set_mu(const double mu) {
      mu_ = mu;
    }

    /** Get mu */
    double get_mu() const {
      return mu_;
    }

    /**
     * Get the offset of the panel into the raw image array
     */
    int2 get_raw_image_offset() const {
      return raw_image_offset_;
    }

    /**
     * Set the offset of the panel into the raw image array
     */
    void set_raw_image_offset(int2 raw_image_offset) {
      raw_image_offset_ = raw_image_offset;
    }

    /** Get the mask array */
    scitbx::af::shared<int4> get_mask() const {
      scitbx::af::shared<int4> result((scitbx::af::reserve(mask_.size())));
      std::copy(mask_.begin(), mask_.end(), std::back_inserter(result));
      return result;
    }

    /** Set the mask */
    void set_mask(const scitbx::af::const_ref<int4> &mask) {
      mask_.clear();
      std::copy(mask.begin(), mask.end(), std::back_inserter(mask_));
    }

    /** Add an element to the mask */
    void add_mask(int f0, int s0, int f1, int s1) {
      mask_.push_back(int4(f0, s0, f1, s1));
    }

    /** @returns True/False this is the same as the other */
    bool operator==(const PanelData &rhs) const {
      return VirtualPanel::operator==(rhs)
             && is_similar_to(rhs, 1e-5, 1e-5, 1e-5, false, false);
    }

    /** @returns True/False this is not the same as the other */
    bool operator!=(const PanelData &rhs) const {
      return !(*this == rhs);
    }

    bool is_similar_to(const PanelData &rhs,
                       double fast_axis_tolerance,
                       double slow_axis_tolerance,
                       double origin_tolerance,
                       bool static_only,
                       bool ignore_trusted_range = false) const {
      bool result =
        image_size_.const_ref().all_eq(rhs.image_size_.const_ref())
        && pixel_size_.const_ref().all_approx_equal(rhs.pixel_size_.const_ref(), 1e-7);
      if (!static_only) {
        result = result
                 && get_fast_axis().const_ref().all_approx_equal(
                   rhs.get_fast_axis().const_ref(), fast_axis_tolerance)
                 && get_slow_axis().const_ref().all_approx_equal(
                   rhs.get_slow_axis().const_ref(), slow_axis_tolerance)
                 && get_origin().const_ref().all_approx_equal(
                   rhs.get_origin().const_ref(), origin_tolerance);
        if (!ignore_trusted_range) {
          result = result
                   && trusted_range_.const_ref().all_approx_equal(
                     rhs.trusted_range_.const_ref(), 1e-7);
        }
      }
      return result;
    }

  protected:
    tiny<double, 2> pixel_size_;
    tiny<std::size_t, 2> image_size_;
    tiny<double, 2> trusted_range_;
    double thickness_;
    std::string material_;
    double mu_;
    int2 raw_image_offset_;
    scitbx::af::shared<int4> mask_;
  };

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_PANEL_H
