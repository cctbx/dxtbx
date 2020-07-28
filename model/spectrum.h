/*
 * spectrum.h
 */
#ifndef DXTBX_MODEL_SPECTRUM_H
#define DXTBX_MODEL_SPECTRUM_H

#include <iostream>
#include <scitbx/array_family/shared.h>
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
        : energies_(energies),
          weights_(weights) {}

    virtual ~Spectrum() {}

    /* Get the spectrum energies (eV) */
    vecd get_energies_eV() const {
      return energies_;
    }

    /* Get the spectrum weights (unitless) */
    vecd get_weights() const {
      return weights_;
    }

    double get_weighted_energy_eV() const {
      if (energies_.size() == 0)
        return 0;
      double weighted_sum = 0;
      double summed_weights = 0;
      for (size_t i = 0; i < energies_.size(); i++) {
        weighted_sum += energies_[i] * weights_[i];
        summed_weights += weights_[i];
      }
      DXTBX_ASSERT(weighted_sum > 0 && summed_weights > 0);
      return weighted_sum / summed_weights;
    }

    double get_weighted_wavelength() const {
      return 12398.4187 / get_weighted_energy_eV(); //  eV per Å conversion factor
    }

    friend std::ostream &operator<<(std::ostream &os, const Spectrum &s);

  private:
    vecd energies_;
    vecd weights_;
  };

  /** Print Spectrum information */
  inline std::ostream &operator<<(std::ostream &os, const Spectrum &s) {
    os << "Spectrum:\n";
    os << "    weighted wavelength: " << s.get_weighted_wavelength() << "\n";
    return os;
  }

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_BEAM_H
