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
#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/vec3.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/sequence.h>
#include <dxtbx/model/crystal.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  /**
   * A class to represent what's in an experiment.
   *
   * Contains:
   *   - imageset Access to the image data 
   *    (ImageSet, ImageSetLazy, ImageSequence, 
   *     ImageGrid, ImageSequence)
   *   - beam The beam model (MonochromaticBeam, TOFBeam)
   *   - detector The detector model
   *   - goniometer The goniometer model
   *   - sequence The sequence model (Sequence, Scan, TOFSequence)
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
    Experiment(boost::python::object beam,
               boost::shared_ptr<Detector> detector,
               boost::shared_ptr<Goniometer> goniometer,
               boost::python::object sequence,
               boost::shared_ptr<CrystalBase> crystal,
               boost::python::object profile,
               boost::python::object imageset,
               boost::python::object scaling_model,
               std::string identifier)
        : beam_(beam),
          detector_(detector),
          goniometer_(goniometer),
          sequence_(sequence),
          crystal_(crystal),
          profile_(profile),
          imageset_(imageset),
          scaling_model_(scaling_model),
          identifier_(identifier) {}

    /**
     * Check if the beam model is the same.
     */
    bool contains_beam(const boost::python::object &beam) const {
      return beam_ == beam;
    }

    /**
     * Check if the detector model is the same.
     */
    bool contains(const boost::shared_ptr<Detector> &detector) const {
      return detector_ == detector;
    }

    /**
     * Check if the detector model is the same.
     */
    bool contains(const boost::shared_ptr<Goniometer> &goniometer) const {
      return goniometer_ == goniometer;
    }

    /**
     * Check if the goniometer model is the same.
     */
    bool contains_sequence(const boost::python::object &sequence) const {
      return sequence_ == sequence;
    }

    /**
     * Check if the crystal model is the same.
     */
    bool contains(const boost::shared_ptr<CrystalBase> &crystal) const {
      return crystal_ == crystal;
    }

    /**
     * Check models are the same.
     */
    bool contains(boost::python::object obj) const {
      boost::python::extract<boost::python::object > get_beam(obj);
      boost::python::extract<boost::shared_ptr<Detector> > get_detector(obj);
      boost::python::extract<boost::shared_ptr<Goniometer> > get_goniometer(obj);
      boost::python::extract<boost::python::object > get_sequence(obj);
      boost::python::extract<boost::shared_ptr<CrystalBase> > get_crystal(obj);
      if (get_beam.check()) {
        return contains_beam(get_beam());
      } else if (get_detector.check()) {
        return contains(get_detector());
      } else if (get_goniometer.check()) {
        return contains(get_goniometer());
      } else if (get_sequence.check()) {
        return contains_sequence(get_sequence());
      } else if (get_crystal.check()) {
        return contains(get_crystal());
      }
      return profile_ == obj || imageset_ == obj || scaling_model_ == obj;
    }

    /**
     * Compare this experiment with another
     */
    bool operator==(const Experiment &other) const {
      return imageset_ == other.imageset_ && beam_ == other.beam_
             && detector_ == other.detector_ && goniometer_ == other.goniometer_
             && sequence_ == other.sequence_ && profile_ == other.profile_
             && scaling_model_ == other.scaling_model_
             && identifier_ == other.identifier_;
    }

    /**
     * Check that the experiment is consistent
     */
    bool is_consistent() const {
      return true;  // FIXME
    }

    /**
     * Check if this experiment represents a still image
     */
    bool is_still() const {
      if (!goniometer_){return true;} 
      if (!sequence_){return true;} 
      std::string sequence_type = boost::python::extract<std::string>(sequence_.attr("__class__").attr("__name__"));
      if (sequence_type == "Scan"){
        Scan scan = boost::python::extract<Scan>(sequence_);
        return scan.get_oscillation()[1] == 0.0;
      }
      return false;
    }

    /**
     * Check if this experiment represents swept rotation image(s)
     */
    bool is_sequence() const {
      return !is_still();
    }

    bool is_tof_experiment() const{
      if (!sequence_){return false;}
      std::string sequence_type = boost::python::extract<std::string>(sequence_.attr("__class__").attr("__name__"));
      if (sequence_type != "TOFSequence"){
        return false;
      }
      std::string beam_type = boost::python::extract<std::string>(beam_.attr("__class__").attr("__name__"));
      if (beam_type != "TOFBeam"){
        return false;
      }
      std::string imageset_type = boost::python::extract<std::string>(imageset_.attr("__class__").attr("__name__"));
      if(imageset_type != "ImageSequence"){
        return false;
      }
      return true;
    }

    /**
     * Set the beam model
     */
    void set_beam(boost::python::object beam) {
      beam_ = beam;
    }

    /**
     * Get the beam model
     */
    boost::python::object get_beam() const {
      return beam_;
    }

    /**
     * Get the detector model
     */
    void set_detector(boost::shared_ptr<Detector> detector) {
      detector_ = detector;
    }

    /**
     * Get the detector model
     */
    boost::shared_ptr<Detector> get_detector() const {
      return detector_;
    }

    /**
     * Get the goniometer model
     */
    void set_goniometer(boost::shared_ptr<Goniometer> goniometer) {
      goniometer_ = goniometer;
    }

    /**
     * Get the goniometer model
     */
    boost::shared_ptr<Goniometer> get_goniometer() const {
      return goniometer_;
    }

    /**
     * Get the sequence model
     */
    void set_sequence(boost::python::object sequence) {
      sequence_ = sequence;
    }

    /**
     * Get the sequence model
     */
    boost::python::object get_sequence() const {
      return sequence_;
    }

    /**
     * Get the crystal model
     */
    void set_crystal(boost::shared_ptr<CrystalBase> crystal) {
      crystal_ = crystal;
    }

    /**
     * Get the crystal model
     */
    boost::shared_ptr<CrystalBase> get_crystal() const {
      return crystal_;
    }

    /**
     * Get the profile model
     */
    void set_profile(boost::python::object profile) {
      profile_ = profile;
    }

    /**
     * Get the profile model
     */
    boost::python::object get_profile() const {
      return profile_;
    }

    /**
     * Get the imageset model
     */
    void set_imageset(boost::python::object imageset) {
      imageset_ = imageset;
    }

    /**
     * Get the imageset model
     */
    boost::python::object get_imageset() const {
      return imageset_;
    }

    /**
     * Set the scaling model
     */
    void set_scaling_model(boost::python::object scaling_model) {
      scaling_model_ = scaling_model;
    }

    /**
     * Get the scaling model
     */
    boost::python::object get_scaling_model() const {
      return scaling_model_;
    }

    /**
     * Set the identifier
     */
    void set_identifier(std::string identifier) {
      identifier_ = identifier;
    }

    /**
     * Get the identifier
     */
    std::string get_identifier() const {
      return identifier_;
    }

  protected:
    boost::python::object beam_;
    boost::shared_ptr<Detector> detector_;
    boost::shared_ptr<Goniometer> goniometer_;
    boost::python::object sequence_;
    boost::shared_ptr<CrystalBase> crystal_;
    boost::python::object profile_;
    boost::python::object imageset_;
    boost::python::object scaling_model_;
    std::string identifier_;
  };

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_EXPERIMENT_H
