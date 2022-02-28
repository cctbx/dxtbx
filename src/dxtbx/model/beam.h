/*
 * beam.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_BEAM_H
#define DXTBX_MODEL_BEAM_H

#include <iostream>
#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include "model_helpers.h"

namespace dxtbx { namespace model {

  using scitbx::vec3;

  /** Base class for beam objects */
  class BeamBase {
  public:
    virtual ~BeamBase() {}

    virtual vec3<double> get_sample_to_source_direction() const = 0;
    virtual double get_wavelength() const = 0;
    virtual double get_divergence() const = 0;
    // Get the standard deviation of the beam divergence
    virtual double get_sigma_divergence() const = 0;
    // Get the wave vector in units of inverse angstroms
    virtual vec3<double> get_s0() const = 0;
    virtual vec3<double> get_unit_s0() const = 0;
    virtual vec3<double> get_polarization_normal() const = 0;
    virtual double get_polarization_fraction() const = 0;
    virtual double get_flux() const = 0;
    virtual double get_transmission() const = 0;
    virtual std::size_t get_num_scan_points() const = 0;
    virtual scitbx::af::shared<vec3<double> > get_s0_at_scan_points() const = 0;
    virtual vec3<double> get_s0_at_scan_point(std::size_t index) const = 0;

    virtual void set_direction(vec3<double> direction) = 0;
    virtual void set_wavelength(double wavelength) = 0;
    // Set the wave vector in units of inverse angstroms
    virtual void set_s0(vec3<double> s0) = 0;
    virtual void set_unit_s0(vec3<double> unit_s0) = 0;
    virtual void set_divergence(double divergence) = 0;
    // Set the standard deviation of the beam divergence
    virtual void set_sigma_divergence(double sigma_divergence) = 0;
    virtual void set_polarization_normal(vec3<double> polarization_normal) = 0;
    virtual void set_polarization_fraction(double polarization_fraction) = 0;
    virtual void set_flux(double flux) = 0;
    virtual void set_transmission(double transmission) = 0;
    virtual void set_s0_at_scan_points(
      const scitbx::af::const_ref<vec3<double> > &s0) = 0;

    virtual void reset_scan_points() = 0;
    virtual bool is_similar_to(const BeamBase &rhs,
                               double wavelength_tolerance,
                               double direction_tolerance,
                               double polarization_normal_tolerance,
                               double polarization_fraction_tolerance) const = 0;
    virtual void rotate_around_origin(vec3<double> axis, double angle) = 0;
    virtual bool operator!=(const BeamBase &rhs) const = 0;
    virtual bool operator==(const BeamBase &rhs) const = 0;
  };

  /** A class to represent a simple beam. */
  class Beam : public BeamBase {
  public:
    Beam()
        : wavelength_(0.0),
          direction_(0.0, 0.0, 1.0),
          divergence_(0.0),
          sigma_divergence_(0.0),
          polarization_normal_(0.0, 1.0, 0.0),
          polarization_fraction_(0.999),
          flux_(0),
          transmission_(1.0) {}

    /**
     * @param s0 The incident beam vector.
     */
    Beam(vec3<double> s0)
        : divergence_(0.0),
          sigma_divergence_(0.0),
          polarization_normal_(0.0, 1.0, 0.0),
          polarization_fraction_(0.999),
          flux_(0),
          transmission_(1.0) {
      DXTBX_ASSERT(s0.length() > 0);
      wavelength_ = 1.0 / s0.length();
      direction_ = -s0.normalize();
    }

    /**
     * @param direction The beam direction vector from sample to source
     * @param wavelength The wavelength of the beam
     */
    Beam(vec3<double> direction, double wavelength)
        : wavelength_(wavelength),
          divergence_(0.0),
          sigma_divergence_(0.0),
          polarization_normal_(0.0, 1.0, 0.0),
          polarization_fraction_(0.999),
          flux_(0),
          transmission_(1.0) {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
    }

    /**
     * @param s0 The incident beam vector.
     * @param divergence The beam divergence
     * @param sigma_divergence The standard deviation of the beam divergence
     */
    Beam(vec3<double> s0, double divergence, double sigma_divergence)
        : divergence_(divergence),
          sigma_divergence_(sigma_divergence),
          polarization_normal_(0.0, 1.0, 0.0),
          polarization_fraction_(0.999),
          flux_(0),
          transmission_(1.0) {
      DXTBX_ASSERT(s0.length() > 0);
      wavelength_ = 1.0 / s0.length();
      direction_ = -s0.normalize();
    }

    /**
     * @param direction The beam direction vector from sample to source
     * @param wavelength The wavelength of the beam
     * @param divergence The beam divergence
     * @param sigma_divergence The standard deviation of the beam divergence
     */
    Beam(vec3<double> direction,
         double wavelength,
         double divergence,
         double sigma_divergence)
        : wavelength_(wavelength),
          divergence_(divergence),
          sigma_divergence_(sigma_divergence),
          polarization_normal_(0.0, 1.0, 0.0),
          polarization_fraction_(0.999),
          flux_(0),
          transmission_(1.0) {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
    }

    /**
     * @param direction The beam direction vector from sample to source
     * @param wavelength The wavelength of the beam
     * @param divergence The beam divergence
     * @param sigma_divergence The standard deviation of the beam divergence
     * @param polarization_normal The polarization plane
     * @param polarization_fraction The polarization fraction
     * @param flux The beam flux
     * @param transmission The beam transmission
     */
    Beam(vec3<double> direction,
         double wavelength,
         double divergence,
         double sigma_divergence,
         vec3<double> polarization_normal,
         double polarization_fraction,
         double flux,
         double transmission)
        : wavelength_(wavelength),
          divergence_(divergence),
          sigma_divergence_(sigma_divergence),
          polarization_normal_(polarization_normal),
          polarization_fraction_(polarization_fraction),
          flux_(flux),
          transmission_(transmission) {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
    }

    virtual ~Beam() {}

    vec3<double> get_sample_to_source_direction() const {
      return direction_;
    }

    virtual double get_wavelength() const {
      return wavelength_;
    }

    double get_divergence() const {
      return divergence_;
    }

    /** Get the standard deviation of the beam divergence */
    double get_sigma_divergence() const {
      return sigma_divergence_;
    }

    /** Set the sample to source direction. */
    void set_direction(vec3<double> direction) {
      DXTBX_ASSERT(direction.length() > 0);
      direction_ = direction.normalize();
    }

    virtual void set_wavelength(double wavelength) {
      wavelength_ = wavelength;
    }

    virtual vec3<double> get_s0() const {
      DXTBX_ASSERT(wavelength_ != 0.0);
      return -direction_ * 1.0 / wavelength_;
    }

    virtual void set_s0(vec3<double> s0) {
      DXTBX_ASSERT(s0.length() > 0);
      direction_ = -s0.normalize();
      wavelength_ = 1.0 / s0.length();
    }

    vec3<double> get_unit_s0() const {
      return -direction_;
    }

    void set_unit_s0(vec3<double> unit_s0) {
      DXTBX_ASSERT(unit_s0.length() > 0);
      direction_ = -(unit_s0.normalize());
    }

    void set_divergence(double divergence) {
      divergence_ = divergence;
    }

    /** Set the standard deviation of the beam divergence */
    void set_sigma_divergence(double sigma_divergence) {
      sigma_divergence_ = sigma_divergence;
    }

    vec3<double> get_polarization_normal() const {
      return polarization_normal_;
    }

    double get_polarization_fraction() const {
      return polarization_fraction_;
    }

    void set_polarization_normal(vec3<double> polarization_normal) {
      polarization_normal_ = polarization_normal;
    }

    void set_polarization_fraction(double polarization_fraction) {
      polarization_fraction_ = polarization_fraction;
    }

    void set_flux(double flux) {
      flux_ = flux;
    }

    void set_transmission(double transmission) {
      transmission_ = transmission;
    }

    double get_flux() const {
      return flux_;
    }

    double get_transmission() const {
      return transmission_;
    }

    std::size_t get_num_scan_points() const {
      return s0_at_scan_points_.size();
    }

    void set_s0_at_scan_points(const scitbx::af::const_ref<vec3<double> > &s0) {
      s0_at_scan_points_ = scitbx::af::shared<vec3<double> >(s0.begin(), s0.end());
    }

    scitbx::af::shared<vec3<double> > get_s0_at_scan_points() const {
      return s0_at_scan_points_;
    }

    vec3<double> get_s0_at_scan_point(std::size_t index) const {
      DXTBX_ASSERT(index < s0_at_scan_points_.size());
      return s0_at_scan_points_[index];
    }

    void reset_scan_points() {
      s0_at_scan_points_.clear();
    }

    virtual bool operator==(const BeamBase &rhs) const {
      double eps = 1.0e-6;

      // scan-varying model checks
      if (get_num_scan_points() > 0) {
        if (get_num_scan_points() != rhs.get_num_scan_points()) {
          return false;
        }
        for (std::size_t j = 0; j < get_num_scan_points(); ++j) {
          vec3<double> this_s0 = get_s0_at_scan_point(j);
          vec3<double> other_s0 = rhs.get_s0_at_scan_point(j);
          double d_s0 = 0.0;
          for (std::size_t i = 0; i < 3; ++i) {
            d_s0 += std::abs(this_s0[i] - other_s0[i]);
          }
          if (d_s0 > eps) {
            return false;
          }
        }
      }

      // static model checks
      return std::abs(angle_safe(direction_, rhs.get_sample_to_source_direction()))
               <= eps
             && std::abs(wavelength_ - rhs.get_wavelength()) <= eps
             && std::abs(divergence_ - rhs.get_divergence()) <= eps
             && std::abs(sigma_divergence_ - rhs.get_sigma_divergence()) <= eps
             && std::abs(
                  angle_safe(polarization_normal_, rhs.get_polarization_normal()))
                  <= eps
             && std::abs(polarization_fraction_ - rhs.get_polarization_fraction())
                  <= eps;
    }

    virtual bool is_similar_to(const BeamBase &rhs,
                       double wavelength_tolerance,
                       double direction_tolerance,
                       double polarization_normal_tolerance,
                       double polarization_fraction_tolerance) const {
      // scan varying model checks
      if (get_num_scan_points() != rhs.get_num_scan_points()) {
        return false;
      }
      for (std::size_t i = 0; i < get_num_scan_points(); ++i) {
        vec3<double> s0_a = get_s0_at_scan_point(i);
        vec3<double> s0_b = rhs.get_s0_at_scan_point(i);

        vec3<double> us0_a = s0_a.normalize();
        vec3<double> us0_b = s0_b.normalize();
        if (std::abs(angle_safe(us0_a, us0_b)) > direction_tolerance) {
          return false;
        }

        double wavelength_a = 1.0 / s0_a.length();
        double wavelength_b = 1.0 / s0_b.length();
        if (std::abs(wavelength_a - wavelength_b) > wavelength_tolerance) {
          return false;
        }
      }

      // static model checks
      return std::abs(angle_safe(direction_, rhs.get_sample_to_source_direction()))
               <= direction_tolerance
             && std::abs(wavelength_ - rhs.get_wavelength()) <= wavelength_tolerance
             && std::abs(
                  angle_safe(polarization_normal_, rhs.get_polarization_normal()))
                  <= polarization_normal_tolerance
             && std::abs(polarization_fraction_ - rhs.get_polarization_fraction())
                  <= polarization_fraction_tolerance;
    }

    bool operator!=(const BeamBase &rhs) const {
      return !(*this == rhs);
    }

    void rotate_around_origin(vec3<double> axis, double angle) {
      const double EPS = 1e-7;
      DXTBX_ASSERT(std::abs(direction_ * polarization_normal_) < EPS);
      direction_ = direction_.rotate_around_origin(axis, angle);
      polarization_normal_ = polarization_normal_.rotate_around_origin(axis, angle);
    }

    friend std::ostream &operator<<(std::ostream &os, const Beam &b);

  protected:
    double wavelength_;
    vec3<double> direction_;
    double divergence_;
    double sigma_divergence_;
    vec3<double> polarization_normal_;
    double polarization_fraction_;
    double flux_;
    double transmission_;
    scitbx::af::shared<vec3<double> > s0_at_scan_points_;
  };

  /** Print Beam information */
  inline std::ostream &operator<<(std::ostream &os, const Beam &b) {
    os << "Beam:\n";
    os << "    wavelength: " << b.get_wavelength() << "\n";
    os << "    sample to source direction : "
       << b.get_sample_to_source_direction().const_ref() << "\n";
    os << "    divergence: " << b.get_divergence() << "\n";
    os << "    sigma divergence: " << b.get_sigma_divergence() << "\n";
    os << "    polarization normal: " << b.get_polarization_normal().const_ref()
       << "\n";
    os << "    polarization fraction: " << b.get_polarization_fraction() << "\n";
    os << "    flux: " << b.get_flux() << "\n";
    os << "    transmission: " << b.get_transmission() << "\n";
    return os;
  }

  /** 
   * A class to represent a time-of-flight beam. 
   * Wavelengths are reflection-specific for time-of-flight experments.
   * To satisfy the BeamBase interface TOFBeam can hold data
   * of a given reflection (reflection_tof_ and reflection_s1_length_)
   */
  class TOFBeam : public Beam{
  public:

    TOFBeam()
      : sample_to_moderator_distance_(0){}

    /**
     * @param direction The beam direction vector from sample to source
     * @param sample_to_moderator_distance (m)
     */
    TOFBeam(vec3<double> direction, double sample_to_moderator_distance){
      set_direction(direction);
      set_sample_to_moderator_distance(sample_to_moderator_distance);
    }

    /**
     * @param direction The beam direction vector from sample to source
     * @param sample_to_moderator_distance (m)
     * @param divergence The beam divergence
     * @param sigma_divergence The standard deviation of the beam divergence
     * @param polarization_normal The polarization plane
     * @param polarization_fraction The polarization fraction
     * @param flux The beam flux
     * @param transmission The beam transmission
     */
    TOFBeam(
      vec3<double> direction,
      double sample_to_moderator_distance,
      double divergence,
      double sigma_divergence,
      vec3<double> polarization_normal,
      double polarization_fraction,
      double flux,
      double transmission
      ){
      set_direction(direction);
      set_sample_to_moderator_distance(sample_to_moderator_distance);
      set_divergence(divergence);
      set_sigma_divergence(sigma_divergence);
      set_polarization_normal(polarization_normal);
      set_polarization_fraction(polarization_fraction);
      set_flux(flux);
      set_transmission(transmission);
      }


    virtual ~TOFBeam(){}

    double get_sample_to_moderator_distance() const{
      return sample_to_moderator_distance_;
    }

    boost::optional<double> get_reflection_tof() const{
      return reflection_tof_;
    }

    boost::optional<double> get_reflection_s1_length() const{
      return reflection_s1_length_;
    }

    double get_wavelength() const{
      DXTBX_ASSERT(reflection_tof_);
      DXTBX_ASSERT(reflection_s1_length_);
      return get_wavelength(*reflection_tof_, *reflection_s1_length_);
    }

    /**
     * @param reflection_tof time-of-flight (s)
     * @param reflection_s1_length length of diffracted beam vector (m)
     * @returns wavelength (A)
     */
    double get_wavelength(double reflection_tof, double reflection_s1_length) const{
      // TODO use scitbx/constants.h instead when #738 has been merged
      double Planck = 6.62607015e-34;
      double m_n = 1.67492749804e-27;
      double L = reflection_s1_length + sample_to_moderator_distance_;
      return ((Planck * reflection_tof) / (m_n * L)) * std::pow(10,10);
    }

    vec3<double> get_s0() const{
      return -direction_ * 1.0 / get_wavelength();
    }

    /**
     * @param reflection_tof time-of-flight (s)
     * @param reflection_s1_length length of diffracted beam vector (m)
     * @returns s0 (A^-1)
     */
    vec3<double> get_s0(double reflection_tof, double reflection_s1_length) const{
      return -direction_ * 1.0 / get_wavelength(
                                          reflection_tof, 
                                          reflection_s1_length
                                          );
    }

    /**
     * @param wavelength (A)
     */
    void set_wavelength(double wavelength){
      DXTBX_ERROR("wavelength is set via set_reflection for TOFBeam");
    } 

    /**
     * @param s0 (A^-1)
     */
    void set_s0(vec3<double> s0){
      DXTBX_ERROR("s0 is set via set_reflection for TOFBeam");
    } 

    /**
     * @param sample_to_moderator_distance (m)
     */
    void set_sample_to_moderator_distance(double sample_to_moderator_distance){
      DXTBX_ASSERT(sample_to_moderator_distance > 0);
      sample_to_moderator_distance_ = sample_to_moderator_distance;
    }

    /**
     * @param tof time-of-flight (s)
     * @param s1_length length of diffracted beam vector (m)
     */
    void set_reflection(double tof, double s1_length){
      DXTBX_ASSERT(tof > 0);
      reflection_tof_ = tof;
      DXTBX_ASSERT(s1_length > 0);
      reflection_s1_length_ = s1_length;
    }

    void clear_reflection(){
      reflection_tof_ = boost::none;
      reflection_s1_length_ = boost::none;
    }

    bool has_reflection() const{
      return (reflection_tof_ && reflection_s1_length_);
    }


    virtual bool operator==(const BeamBase &rhs) const{
      double eps = 1.0e-6;
      return std::abs(angle_safe(direction_, rhs.get_sample_to_source_direction()))
               <= eps
             && std::abs(divergence_ - rhs.get_divergence()) <= eps
             && std::abs(sigma_divergence_ - rhs.get_sigma_divergence()) <= eps
             && std::abs(
                  angle_safe(polarization_normal_, rhs.get_polarization_normal()))
                  <= eps
             && std::abs(polarization_fraction_ - rhs.get_polarization_fraction())
                  <= eps;

    }

    virtual bool operator==(const TOFBeam &rhs) const {
      double eps = 1.0e-6;

      // Only compare scan varying and wavelength if both have reflection data
      if (has_reflection() && rhs.has_reflection()){

        if (std::abs(get_wavelength() - rhs.get_wavelength()) > eps){
          return false;
        }

        // scan-varying model checks
        if (get_num_scan_points() > 0) {
          if (get_num_scan_points() != rhs.get_num_scan_points()) {
            return false;
          }
          for (std::size_t j = 0; j < get_num_scan_points(); ++j) {
            vec3<double> this_s0 = get_s0_at_scan_point(j);
            vec3<double> other_s0 = rhs.get_s0_at_scan_point(j);
            double d_s0 = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
              d_s0 += std::abs(this_s0[i] - other_s0[i]);
            }
            if (d_s0 > eps) {
              return false;
            }
          }
        }
      }

      // static model checks
      return std::abs(angle_safe(direction_, rhs.get_sample_to_source_direction()))
               <= eps
             && std::abs(divergence_ - rhs.get_divergence()) <= eps
             && std::abs(sigma_divergence_ - rhs.get_sigma_divergence()) <= eps
             && std::abs(
                  angle_safe(polarization_normal_, rhs.get_polarization_normal()))
                  <= eps
             && std::abs(polarization_fraction_ - rhs.get_polarization_fraction())
                  <= eps;
    }

    bool operator!=(const TOFBeam &rhs) const {
      return !(*this == rhs);
    }

    virtual bool is_similar_to(const BeamBase &rhs,
                       double wavelength_tolerance,
                       double direction_tolerance,
                       double polarization_normal_tolerance,
                       double polarization_fraction_tolerance) const {

      // static model checks
      return std::abs(angle_safe(direction_, rhs.get_sample_to_source_direction()))
               <= direction_tolerance
             && std::abs(
                  angle_safe(polarization_normal_, rhs.get_polarization_normal()))
                  <= polarization_normal_tolerance
             && std::abs(polarization_fraction_ - rhs.get_polarization_fraction())
                  <= polarization_fraction_tolerance;

    }

    virtual bool is_similar_to(const TOFBeam &rhs,
                       double wavelength_tolerance,
                       double direction_tolerance,
                       double polarization_normal_tolerance,
                       double polarization_fraction_tolerance) const {


      // Only compare scan varying and wavelength if both have reflection data
      if (has_reflection() && rhs.has_reflection()){

        if (std::abs(get_wavelength() - rhs.get_wavelength()) > wavelength_tolerance){
          return false;
        }

        // scan varying model checks
        if (get_num_scan_points() != rhs.get_num_scan_points()) {
          return false;
        }
        for (std::size_t i = 0; i < get_num_scan_points(); ++i) {
          vec3<double> s0_a = get_s0_at_scan_point(i);
          vec3<double> s0_b = rhs.get_s0_at_scan_point(i);

          vec3<double> us0_a = s0_a.normalize();
          vec3<double> us0_b = s0_b.normalize();
          if (std::abs(angle_safe(us0_a, us0_b)) > direction_tolerance) {
            return false;
          }

          double wavelength_a = 1.0 / s0_a.length();
          double wavelength_b = 1.0 / s0_b.length();
          if (std::abs(wavelength_a - wavelength_b) > wavelength_tolerance) {
            return false;
          }
        }
      }

      // static model checks
      return std::abs(angle_safe(direction_, rhs.get_sample_to_source_direction()))
               <= direction_tolerance
             && std::abs(
                  angle_safe(polarization_normal_, rhs.get_polarization_normal()))
                  <= polarization_normal_tolerance
             && std::abs(polarization_fraction_ - rhs.get_polarization_fraction())
                  <= polarization_fraction_tolerance;
    }

    friend std::ostream &operator<<(std::ostream &os, const TOFBeam &b);


  private:
    double sample_to_moderator_distance_;
    boost::optional<double> reflection_tof_;
    boost::optional<double> reflection_s1_length_;
  };

  /** Print TOFBeam information */
  inline std::ostream &operator<<(std::ostream &os, const TOFBeam &b) {
    os << "TOFBeam:\n";
    os << "    sample to source direction : "
       << b.get_sample_to_source_direction().const_ref() << "\n";
    os << "    sample_to_moderator_distance: " << b.get_sample_to_moderator_distance() << "\n";
    os << "    divergence: " << b.get_divergence() << "\n";
    os << "    sigma divergence: " << b.get_sigma_divergence() << "\n";
    os << "    polarization normal: " << b.get_polarization_normal().const_ref()
       << "\n";
    os << "    polarization fraction: " << b.get_polarization_fraction() << "\n";
    os << "    flux: " << b.get_flux() << "\n";
    os << "    transmission: " << b.get_transmission() << "\n";
    return os;
  }

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_BEAM_H
