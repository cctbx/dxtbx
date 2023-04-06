/*
 * scan_helpers.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_SCAN_HELPERS_H
#define DXTBX_MODEL_SCAN_HELPERS_H

#include <cmath>
#include <limits>
#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/array_family/shared.h>
#include <dxtbx/error.h>
#include <iostream>

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::constants::two_pi;
  using std::floor;

  /** Convert the angle mod 2PI */
  inline double mod_2pi(double angle) {
    return angle - two_pi * floor(angle / two_pi);
  }

  /**
   * Check if the angle is within the given range. The angular
   * range can be any two angles, plus or minus. The angle is to check can
   * also be any angle. The angle is considered within the range if the range
   * spans more than 2PI degrees or the angle is within the two range angles
   * when mod 2PI.
   * @param range The angular range
   * @param angle The angle to check
   * @returns True/False the angle is within the range
   */

  inline bool is_angle_in_range(vec2<double> range, double angle) {
    auto sign_changed = [](double old_value, double new_value) {
      if (old_value < 0 && new_value > 0) {
        return true;
      }
      if (old_value > 0 && new_value < 0) {
        return true;
      }
      return false;
    };

    // All angles in range mod_2pi
    if (std::abs(range[1] - range[0]) >= two_pi) {
      return true;
    }

    double range_start;
    double range_end;
    double eps = std::numeric_limits<double>::epsilon();

    // Account for range possibly being negative
    range_start = std::min(range[0], range[1]);
    range_end = std::max(range[0], range[1]);

    // Simple case
    if (angle >= range_start and angle <= range_end) {
      return true;
    }

    // Case where angle and/or range is outside of 0,360
    float angle_mod_2pi = mod_2pi(angle);
    float range0_mod_2pi = mod_2pi(range[0]);
    float range1_mod_2pi = mod_2pi(range[1]);

    // If either range element is 0, and the other element has changed sign, set to 2pi
    if (range[0] <= eps && range[0] >= -eps && sign_changed(range[1], range1_mod_2pi)) {
      range0_mod_2pi = two_pi;
    }
    if (range[1] <= eps && range[1] >= -eps && sign_changed(range[0], range0_mod_2pi)) {
      range1_mod_2pi = two_pi;
    }

    // Account for range originally being negative
    range_start = std::min(range0_mod_2pi, range1_mod_2pi);
    range_end = std::max(range0_mod_2pi, range1_mod_2pi);
    if (angle_mod_2pi >= range_start and angle_mod_2pi <= range_end) {
      return true;
    }
    return false;
  }

  /**
   * A function to get the range of equivalent angles mod 2PI degrees that lie
   * in the given angular range which can be more than 2PI degrees. For example
   * if the range is given as (A, B), the returned value will be (a, b) where
   * A <= a < A + 2PI; B - 2PI < b < B.
   * @param angle The angle to check
   * @returns A pair of angles (a, b) where a = angle + n2PI and
   * b = angle + m2pi and both lie within the caches range values. If
   * b < a, the range is invalid.
   */
  inline vec2<double> get_range_of_mod2pi_angles(vec2<double> range, double angle) {
    return vec2<double>(angle - two_pi * floor((angle - range[0]) / two_pi),
                        angle + two_pi * floor((range[1] - angle) / two_pi));
  }

  /**
   * A function to get the all the angles mod 2pi a given angle that lie
   * in the given angular range which can be more than 2pi degrees.
   * Calculate and return an array of angles. If no angles are in the range,
   * then the array is empty.
   * @param angle The angle to use.
   * @returns An array of angles, a, where a = a + n2pi.
   */
  inline scitbx::af::shared<double> get_mod2pi_angles_in_range(vec2<double> range,
                                                               double angle) {
    scitbx::af::shared<double> result;
    vec2<double> angle_range = get_range_of_mod2pi_angles(range, angle);
    int n_angles = 1 + (int)floor((angle_range[1] - angle_range[0]) / two_pi);
    if (n_angles > 0) {
      result.resize(n_angles);
      for (std::size_t i = 0; i < n_angles; ++i) {
        result[i] = angle_range[0] + i * two_pi;
      }
    }
    return result;
  }

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_SCAN_HELPERS_H
