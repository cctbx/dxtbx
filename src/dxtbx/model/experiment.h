/*
 * experiment.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_EXPERIMENT_H
#define DXTBX_MODEL_EXPERIMENT_H

#include <iostream>
#include <cmath>
#include <memory>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/vec3.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  enum ExperimentType { rotation = 1, still = 2, tof = 3 };

  /**
   * A class to represent what's in an experiment.
   *
   * Contains:
   *   - imageset Access to the image data
   *   - beam The beam model
   *   - detector The detector model
   *   - goniometer The goniometer model
   *   - scan The scan model
   *   - crystal The crystal model
   *   - profile The profile model
   *   - scaling_model The scaling model
   *
   * Some of these may be set to "None"
   *
   */
  class Experiment {
  public:
    Experiment() {}

    /**
     * Initialise the experiment with models
     */
    Experiment(std::shared_ptr<BeamBase> beam,
               std::shared_ptr<Detector> detector,
               std::shared_ptr<Goniometer> goniometer,
               std::shared_ptr<Scan> scan,
               std::shared_ptr<CrystalBase> crystal,
               boost::python::object profile,
               boost::python::object imageset,
               boost::python::object scaling_model,
               std::string identifier)
        : beam_(beam),
          detector_(detector),
          goniometer_(goniometer),
          scan_(scan),
          crystal_(crystal),
          profile_(profile),
          imageset_(imageset),
          scaling_model_(scaling_model),
          identifier_(identifier) {}

    /**
     * Check if the beam model is the same.
     */
    bool contains(const std::shared_ptr<BeamBase> &beam) const {
      return beam_ == beam;
    }

    /**
     * Check if the detector model is the same.
     */
    bool contains(const std::shared_ptr<Detector> &detector) const {
      return detector_ == detector;
    }

    /**
     * Check if the detector model is the same.
     */
    bool contains(const std::shared_ptr<Goniometer> &goniometer) const {
      return goniometer_ == goniometer;
    }

    /**
     * Check if the goniometer model is the same.
     */
    bool contains(const std::shared_ptr<Scan> &scan) const {
      return scan_ == scan;
    }

    /**
     * Check if the crystal model is the same.
     */
    bool contains(const std::shared_ptr<CrystalBase> &crystal) const {
      return crystal_ == crystal;
    }

    /**
     * Check models are the same.
     */
    bool contains(boost::python::object obj) const {
      boost::python::extract<std::shared_ptr<BeamBase> > get_beam(obj);
      boost::python::extract<std::shared_ptr<Detector> > get_detector(obj);
      boost::python::extract<std::shared_ptr<Goniometer> > get_goniometer(obj);
      boost::python::extract<std::shared_ptr<Scan> > get_scan(obj);
      boost::python::extract<std::shared_ptr<CrystalBase> > get_crystal(obj);
      if (get_beam.check()) {
        return contains(get_beam());
      } else if (get_detector.check()) {
        return contains(get_detector());
      } else if (get_goniometer.check()) {
        return contains(get_goniometer());
      } else if (get_scan.check()) {
        return contains(get_scan());
      } else if (get_crystal.check()) {
        return contains(get_crystal());
      }
      return profile_ == obj || imageset_ == obj || scaling_model_ == obj;
    }

    bool operator==(const Experiment &other) const {
      return imageset_ == other.imageset_ && beam_ == other.beam_
             && detector_ == other.detector_ && goniometer_ == other.goniometer_
             && scan_ == other.scan_ && profile_ == other.profile_
             && scaling_model_ == other.scaling_model_
             && identifier_ == other.identifier_;
    }

    /**
     * Check if this experiment represents a still image
     */
    [[deprecated("use get_type()==ExperimentType.still instead")]] bool is_still()
      const {
      return !goniometer_ || !scan_ || scan_->is_still();
    }

    /**
     * Check if this experiment represents swept rotation image(s)
     */
    [[deprecated("use get_type()==ExperimentType.rotation instead")]] bool is_sequence()
      const {
      return !is_still();
    }

    ExperimentType get_type() {
      if (scan_ && scan_->contains("time_of_flight")) {
        return tof;
      }
      if (!goniometer_ || !scan_ || scan_->is_still()) {
        return still;
      } else {
        return rotation;
      }
    }

    bool is_consistent() const {
      return true;  // FIXME
    }

    void set_beam(std::shared_ptr<BeamBase> beam) {
      beam_ = beam;
    }

    std::shared_ptr<BeamBase> get_beam() const {
      return beam_;
    }

    void set_detector(std::shared_ptr<Detector> detector) {
      detector_ = detector;
    }

    std::shared_ptr<Detector> get_detector() const {
      return detector_;
    }

    void set_goniometer(std::shared_ptr<Goniometer> goniometer) {
      goniometer_ = goniometer;
    }

    std::shared_ptr<Goniometer> get_goniometer() const {
      return goniometer_;
    }

    void set_scan(std::shared_ptr<Scan> scan) {
      scan_ = scan;
    }

    std::shared_ptr<Scan> get_scan() const {
      return scan_;
    }

    void set_crystal(std::shared_ptr<CrystalBase> crystal) {
      crystal_ = crystal;
    }

    std::shared_ptr<CrystalBase> get_crystal() const {
      return crystal_;
    }

    void set_profile(boost::python::object profile) {
      profile_ = profile;
    }

    boost::python::object get_profile() const {
      return profile_;
    }

    void set_imageset(boost::python::object imageset) {
      imageset_ = imageset;
    }

    boost::python::object get_imageset() const {
      return imageset_;
    }

    void set_scaling_model(boost::python::object scaling_model) {
      scaling_model_ = scaling_model;
    }

    boost::python::object get_scaling_model() const {
      return scaling_model_;
    }

    void set_identifier(std::string identifier) {
      identifier_ = identifier;
    }

    std::string get_identifier() const {
      return identifier_;
    }

  protected:
    std::shared_ptr<BeamBase> beam_;
    std::shared_ptr<Detector> detector_;
    std::shared_ptr<Goniometer> goniometer_;
    std::shared_ptr<Scan> scan_;
    std::shared_ptr<CrystalBase> crystal_;
    boost::python::object profile_;
    boost::python::object imageset_;
    boost::python::object scaling_model_;
    std::string identifier_;
  };

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_EXPERIMENT_H
