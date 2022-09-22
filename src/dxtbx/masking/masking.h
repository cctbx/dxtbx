/*
 * masking.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MASKING_H
#define DXTBX_MASKING_H

#include <algorithm>
#include <scitbx/array_family/shared.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/panel.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace masking {

  using dxtbx::model::BeamBase;
  using dxtbx::model::Panel;
  using scitbx::vec2;
  using scitbx::vec3;

  /**
   * Function to mask rectangle
   * @param mask The mask array
   * @param x0 The low x range
   * @param x1 The high x range
   * @param y0 The low y range
   * @param y1 The high y range
   */
  void mask_untrusted_rectangle(scitbx::af::ref<bool, scitbx::af::c_grid<2> > mask,
                                std::size_t x0,
                                std::size_t x1,
                                std::size_t y0,
                                std::size_t y1) {
    std::size_t height = mask.accessor()[0];
    std::size_t width = mask.accessor()[1];
    DXTBX_ASSERT(x1 > x0);
    DXTBX_ASSERT(y1 > y0);
    DXTBX_ASSERT(x1 <= width);
    DXTBX_ASSERT(y1 <= height);
    for (std::size_t j = y0; j < y1; ++j) {
      for (std::size_t i = x0; i < x1; ++i) {
        mask(j, i) = false;
      }
    }
  }

  /**
   * Function to mask a circle
   * @param mask The mask array
   * @param x0 The low x range
   * @param x1 The high x range
   * @param y0 The low y range
   * @param y1 The high y range
   */
  void mask_untrusted_circle(scitbx::af::ref<bool, scitbx::af::c_grid<2> > mask,
                             double xc,
                             double yc,
                             double radius) {
    DXTBX_ASSERT(radius > 0);
    std::size_t height = mask.accessor()[0];
    std::size_t width = mask.accessor()[1];
    int x0 = (int)std::floor(xc - radius);
    int y0 = (int)std::floor(yc - radius);
    int x1 = (int)std::ceil(xc + radius);
    int y1 = (int)std::ceil(yc + radius);
    x0 = std::max(x0, 0);
    y0 = std::max(y0, 0);
    x1 = std::min(x1, (int)width);
    y1 = std::min(y1, (int)height);
    DXTBX_ASSERT(x1 > x0);
    DXTBX_ASSERT(y1 > y0);
    double r2 = radius * radius;
    for (std::size_t j = y0; j < y1; ++j) {
      for (std::size_t i = x0; i < x1; ++i) {
        if ((i - xc) * (i - xc) + (j - yc) * (j - yc) < r2) {
          mask(j, i) = false;
        }
      }
    }
  }

  /**
   * Check point is inside polygon
   * @param poly - The polygon
   * @param x The x coord
   * @param y The y coord
   * @returns True/False the point is inside the polygon
   */
  bool is_inside_polygon(const scitbx::af::const_ref<vec2<double> > &poly,
                         double x,
                         double y) {
    // http://en.wikipedia.org/wiki/Point_in_polygon
    // http://en.wikipedia.org/wiki/Even-odd_rule
    std::size_t num = poly.size();
    std::size_t j = num - 1;
    bool inside = false;
    for (std::size_t i = 0; i < num; ++i) {
      if (((poly[i][1] > y) != (poly[j][1] > y))
          && (x
              < (poly[j][0] - poly[i][0]) * (y - poly[i][1]) / (poly[j][1] - poly[i][1])
                  + poly[i][0])) {
        inside = !inside;
      }
      j = i;
    }
    return inside;
  }

  /**
   * Apply a polygon mask
   * @param mask The mask array
   * @param polygon The polygon
   */
  void mask_untrusted_polygon(scitbx::af::ref<bool, scitbx::af::c_grid<2> > mask,
                              const scitbx::af::const_ref<vec2<double> > &polygon) {
    DXTBX_ASSERT(polygon.size() > 3);
    std::size_t height = mask.accessor()[0];
    std::size_t width = mask.accessor()[1];
    int x0 = (int)std::floor(polygon[0][0]);
    int y0 = (int)std::floor(polygon[0][1]);
    int x1 = x0;
    int y1 = y0;
    for (std::size_t i = 1; i < polygon.size(); ++i) {
      int x = (int)std::floor(polygon[i][0]);
      int y = (int)std::floor(polygon[i][1]);
      if (x < x0) x0 = x;
      if (y < y0) y0 = y;
      if (x > x1) x1 = x;
      if (y > y1) y1 = y;
    }
    x0 = std::max(x0, 0);
    y0 = std::max(y0, 0);
    x1 = std::min(x1 + 1, (int)width);
    y1 = std::min(y1 + 1, (int)height);
    DXTBX_ASSERT(x0 < x1);
    DXTBX_ASSERT(y0 < y1);
    for (std::size_t j = y0; j < y1; ++j) {
      for (std::size_t i = x0; i < x1; ++i) {
        if (is_inside_polygon(polygon, i + 0.5, j + 0.5)) {
          mask(j, i) = false;
        }
      }
    }
  }

  /**
   * Function to add a resolution range mask
   * @param mask The mask array
   * @param beam The beam model
   * @param panel The panel model
   * @param d_min The high resolution limit
   * @param d_max The low resolution limit
   */
  void mask_untrusted_resolution_range(
    scitbx::af::ref<bool, scitbx::af::c_grid<2> > mask,
    const BeamBase &beam,
    const Panel &panel,
    double d_min,
    double d_max) {
    DXTBX_ASSERT(d_min < d_max);
    std::size_t width = panel.get_image_size()[0];
    std::size_t height = panel.get_image_size()[1];
    DXTBX_ASSERT(height == mask.accessor()[0]);
    DXTBX_ASSERT(width == mask.accessor()[1]);
    vec3<double> s0 = beam.get_s0();
    for (std::size_t j = 0; j < height; ++j) {
      for (std::size_t i = 0; i < width; ++i) {
        vec2<double> px(i + 0.5, j + 0.5);
        double d = panel.get_resolution_at_pixel(s0, px);
        if (d_min <= d && d <= d_max) {
          mask(j, i) = false;
        }
      }
    }
  }

  /**
   * A class to mask multiple resolution ranges
   */
  class ResolutionMaskGenerator {
  public:
    /**
     * Initialise the resolution at each pixel
     * @param beam The beam model
     * @param panel The panel model
     */
    ResolutionMaskGenerator(const BeamBase &beam, const Panel &panel)
        : resolution_(
          scitbx::af::c_grid<2>(panel.get_image_size()[1], panel.get_image_size()[0])) {
      vec3<double> s0 = beam.get_s0();
      for (std::size_t j = 0; j < resolution_.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < resolution_.accessor()[1]; ++i) {
          vec2<double> px(i + 0.5, j + 0.5);
          try {
            resolution_(j, i) = panel.get_resolution_at_pixel(s0, px);
          } catch (dxtbx::error const &) {
            // Known failure: resolution at beam center is undefined
            resolution_(j, i) = 0.0;
          }
        }
      }
    }

    /**
     * Apply the mask
     * @param mask The mask
     * @param d_min The high resolution of the range
     * @param d_max The low resolution of the range
     */
    void apply(scitbx::af::ref<bool, scitbx::af::c_grid<2> > mask,
               double d_min,
               double d_max) const {
      DXTBX_ASSERT(d_min < d_max);
      DXTBX_ASSERT(resolution_.accessor()[0] == mask.accessor()[0]);
      DXTBX_ASSERT(resolution_.accessor()[1] == mask.accessor()[1]);
      for (std::size_t j = 0; j < resolution_.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < resolution_.accessor()[1]; ++i) {
          double d = resolution_(j, i);
          if (d_min <= d && d <= d_max) {
            mask(j, i) = false;
          }
        }
      }
    }

  private:
    scitbx::af::versa<double, scitbx::af::c_grid<2> > resolution_;
  };

}}  // namespace dxtbx::masking

#endif /* DXTBX_MASKING_H */
