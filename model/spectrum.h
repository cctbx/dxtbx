/*
 * spectrum.h
 */
#ifndef DXTBX_MODEL_SPECTRUM_H
#define DXTBX_MODEL_SPECTRUM_H

#include <iostream>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <scitbx/constants.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {
  typedef scitbx::af::shared<double> vecd;

  /** A class to represent a 2D spectrum. */
  class Spectrum {
  public:
    /** Default constructor: initialise all to zero */
    Spectrum() {}

    /**
     * Initialise all the spectrum parameters.
     * @param energies The spectrum energies (eV)
     * @param energies The spectrum weights (unitless)
     */
    Spectrum(vecd energies, vecd weights)
        : energies_(energies), weights_(weights), emin_(0.0), emax_(0.0) {
      compute_weighted_energy();
    }

    virtual ~Spectrum() {}

    /* Get the spectrum energies (eV) */
    vecd get_energies_eV() const {
      return energies_;
    }

    /* Get the spectrum weights (unitless) */
    vecd get_weights() const {
      return weights_;
    }

    /* Helper function to compute bandwidth - assumes that the input spectrum
       is somewhat background subtracted */

    void bandwidth_98_percent() {
      if (energies_.size() == 0) return;

      std::vector<double> cdf;
      double total = 0;

      /* Compute CDF */
      for (size_t i = 0; i < energies_.size(); i++) {
        total += weights_[i];
        cdf.push_back(total);
      }

      /* Scan CDF to find 1%, 99% points */
      for (size_t i = 0; i < energies_.size(); i++) {
        if (cdf[i] < 0.01 * total) emin_ = energies_[i];
        if (cdf[i] > 0.99 * total) {
          emax_ = energies_[i];
          break;
        }
      }
    }

    /* Get the bandwidth range */
    double get_emin_eV() {
      if ((emin_ == 0) && (emax_ == 0)) bandwidth_98_percent();
      return emin_;
    }

    /* Get the bandwidth range */
    double get_emax_eV() {
      if ((emin_ == 0) && (emax_ == 0)) bandwidth_98_percent();
      return emax_;
    }

    double get_weighted_energy_eV() const {
      return weighted_energy_;
    }

    double get_weighted_energy_variance() const {
      return weighted_energy_variance_;
    }

    void compute_weighted_energy() {
      if (energies_.size() == 0) {
        weighted_energy_ = 0;
        return;
      }
      double weighted_sum = 0;
      double weighted_sum_sq = 0;
      double summed_weights = 0;
      for (size_t i = 0; i < energies_.size(); i++) {
        weighted_sum += energies_[i] * weights_[i];
        weighted_sum_sq += energies_[i] * energies_[i] * weights_[i];
        summed_weights += weights_[i];
      }
      DXTBX_ASSERT(weighted_sum > 0 && summed_weights > 0);
      weighted_energy_ = weighted_sum / summed_weights;
      weighted_energy_variance_ =
        weighted_sum_sq / summed_weights - (weighted_energy_ * weighted_energy_);
    }

    double get_weighted_wavelength() const {
      return scitbx::constants::factor_ev_angstrom
             / get_weighted_energy_eV();  //  eV per Å conversion factor
    }

    friend std::ostream &operator<<(std::ostream &os, const Spectrum &s);

  private:
    vecd energies_;
    vecd weights_;

    double emin_, emax_;

    double weighted_energy_;
    double weighted_energy_variance_;
  };

  /** Print Spectrum information */
  inline std::ostream &operator<<(std::ostream &os, const Spectrum &s) {
    os << "Spectrum:\n";
    os << "    weighted wavelength: " << s.get_weighted_wavelength() << "\n";
    return os;
  }

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_BEAM_H
