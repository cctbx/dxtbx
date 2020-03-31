/*
 * crystal.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_CRYSTAL_H
#define DXTBX_MODEL_CRYSTAL_H

#include <iostream>
#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <scitbx/constants.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/matrix/move.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/math/angle_derivative.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <cctbx/uctbx.h>
#include <cctbx/crystal_orientation.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  using scitbx::deg_as_rad;
  using scitbx::mat3;
  using scitbx::rad_as_deg;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::double6;

  namespace detail {

    /**
     * Helper function to check if matrix is a rotation matrix
     */
    inline bool is_r3_rotation_matrix(mat3<double> R, double rms_tolerance = 1e-8) {
      mat3<double> rtr = R.transpose() * R;
      mat3<double> identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
      mat3<double> rtrmi = (rtr - identity);
      double norm_sq = 0.0;
      for (std::size_t i = 0; i < 9; ++i) {
        norm_sq += rtrmi[i] * rtrmi[i];
      }
      return std::sqrt(norm_sq) < rms_tolerance
             && std::abs(1.0 - R.determinant()) < rms_tolerance;
    }

    /**
     * Do matrix multiply A * B * A.transpose()
     *
     * ABAT[1] = A[ar, ac] * B[ac, bc] * A[ar, ac].transpose()
     */
    template <typename FloatType>
    FloatType multiply_A_B_AT(const FloatType *A,
                              const FloatType *B,
                              std::size_t ar,
                              std::size_t ac,
                              std::size_t bc) {
      FloatType ABAT = 0;
      FloatType *AB = new FloatType[ar * bc];
      scitbx::matrix::multiply<FloatType, FloatType, FloatType>(A, B, ar, ac, bc, AB);
      scitbx::matrix::multiply_transpose<FloatType, FloatType, FloatType>(
        AB, A, ar, bc, ar, &ABAT);
      delete[] AB;
      return ABAT;
    }

    /**
     * Implement analytical formula of Lefebvre et al. (1999)
     * http://arxiv.org/abs/hep-ex/9909031 to calculate the covariances of elements
     * of mat^-1, given the covariances of mat itself. The input covariance matrix,
     * and the return value of the function, have elements ordered by treating mat as
     * a 1D vector using row-major ordering.
     */
    inline scitbx::af::versa<double, scitbx::af::c_grid<2> >
    matrix_inverse_error_propagation(
      const scitbx::af::ref<double, scitbx::af::c_grid<2> > &mat,
      const scitbx::af::ref<double, scitbx::af::c_grid<2> > &cov_mat) {
      // initialise covariance matrix
      DXTBX_ASSERT(mat.accessor()[0] == mat.accessor()[1]);
      DXTBX_ASSERT(cov_mat.accessor()[0] == cov_mat.accessor()[1]);
      std::size_t n = mat.accessor()[0];
      DXTBX_ASSERT(cov_mat.accessor()[0] == n * n);

      // use flex for nice 2D indexing
      scitbx::af::versa<double, scitbx::af::c_grid<2> > inv_mat(mat.accessor());
      std::copy(mat.begin(), mat.end(), inv_mat.begin());
      scitbx::af::matrix_inversion_in_place(inv_mat.ref());

      scitbx::af::versa<double, scitbx::af::c_grid<2> > inv_cov_mat(cov_mat.accessor(),
                                                                    0.0);
      for (std::size_t alpha = 0; alpha < n; ++alpha) {
        for (std::size_t beta = 0; beta < n; ++beta) {
          for (std::size_t a = 0; a < n; ++a) {
            for (std::size_t b = 0; b < n; ++b) {
              // index into inv_cov_mat after flattening inv_mat
              std::size_t u = alpha * n + beta;
              std::size_t v = a * n + b;

              // skip elements in the lower triangle
              if (v < u) {
                continue;
              }

              // The element u,v of the result is the calculation
              // cov(m^-1[alpha, beta], m^-1[a, b])
              double elt = 0.0;
              for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                  for (std::size_t k = 0; k < n; ++k) {
                    for (std::size_t l = 0; l < n; ++l) {
                      // index into cov_mat after flattening mat
                      std::size_t x = i * n + j;
                      std::size_t y = k * n + l;
                      elt += inv_mat(alpha, i) * inv_mat(j, beta) * inv_mat(a, k)
                             * inv_mat(l, b) * cov_mat(x, y);
                    }
                  }
                }
              }
              inv_cov_mat(u, v) = elt;
            }
          }
        }
      }

      scitbx::matrix::copy_upper_to_lower_triangle_in_place(inv_cov_mat.ref());
      return inv_cov_mat;
    }

  }  // namespace detail

  /** Base class for crystal objects */
  class CrystalBase {
  public:
    virtual ~CrystalBase() {}

    // Set the unit cell parameters
    virtual void set_unit_cell(const vec3<double> &real_space_a,
                               const vec3<double> &real_space_b,
                               const vec3<double> &real_space_c) = 0;
    // Set the unit cell parameters
    virtual void set_unit_cell(const cctbx::uctbx::unit_cell &unit_cell) = 0;
    // Update the B matrix
    virtual void update_B() = 0;
    // Set the U matrix
    virtual void set_U(const mat3<double> &U) = 0;
    // Set the B matrix
    virtual void set_B(const mat3<double> &B) = 0;
    // Set the A matrix
    virtual void set_A(const mat3<double> &A) = 0;
    // @returns the U matrix
    virtual mat3<double> get_U() const = 0;
    // @returns the B matrix
    virtual mat3<double> get_B() const = 0;
    // @returns the A matrix
    virtual mat3<double> get_A() const = 0;
    // @returns the unit_cell
    virtual cctbx::uctbx::unit_cell get_unit_cell() const = 0;
    // @returns the real space vectors
    virtual scitbx::af::shared<vec3<double> > get_real_space_vectors() const = 0;
    // Set the space group
    virtual void set_space_group(cctbx::sgtbx::space_group space_group) = 0;
    // Get the space group
    virtual cctbx::sgtbx::space_group get_space_group() const = 0;
    // @returns the number of scan points
    virtual std::size_t get_num_scan_points() const = 0;
    // Set the A matrix at scan points
    virtual void set_A_at_scan_points(
      const scitbx::af::const_ref<mat3<double> > &A) = 0;
    // Get the A matrix at scan points
    virtual scitbx::af::shared<mat3<double> > get_A_at_scan_points() const = 0;
    //  Get the A matrix at the scan point
    virtual mat3<double> get_A_at_scan_point(std::size_t index) const = 0;
    // Get the U matrix at the scan point
    virtual mat3<double> get_U_at_scan_point(std::size_t index) const = 0;
    // Get the B matrix at the scan point
    virtual mat3<double> get_B_at_scan_point(std::size_t index) const = 0;
    // Get the unit cell at the scan point
    virtual cctbx::uctbx::unit_cell get_unit_cell_at_scan_point(
      std::size_t index) const = 0;
    // Reset the scan points
    virtual void reset_scan_points() = 0;
    // Returns a copy of the current crystal model transformed by the given
    // change of basis operator to the new basis.
    virtual boost::shared_ptr<CrystalBase> change_basis(
      cctbx::sgtbx::change_of_basis_op change_of_basis_op) const = 0;
    // Update crystal with parameters from another model
    virtual void update(const CrystalBase &other) = 0;
    // Rotate the model around an axis and angle
    virtual void rotate_around_origin(vec3<double> axis,
                                      double angle,
                                      bool deg = false) = 0;
    // Check if the models are equal
    virtual bool operator==(const CrystalBase &other) const = 0;
    // Check if models are similar
    virtual bool is_similar_to(const CrystalBase &other,
                               double angle_tolerance = 0.01,
                               double uc_rel_length_tolerance = 0.01,
                               double uc_abs_angle_tolerance = 1.0)
      const = 0;  // FIXME remove these tolerances and move them to other classes?
    // Check if the models are not equal
    virtual bool operator!=(const CrystalBase &other) const = 0;
    // Return the 9*9 covariance matrix of elements of B, if available, otherwise None.
    virtual scitbx::af::versa<double, scitbx::af::c_grid<2> > get_B_covariance()
      const = 0;
    // Set the covariance matrix
    virtual void set_B_covariance(
      const scitbx::af::const_ref<double, scitbx::af::c_grid<2> > &cov) = 0;
    virtual void set_B_covariance_at_scan_points(
      const scitbx::af::const_ref<double, scitbx::af::c_grid<3> > &cov) = 0;
    virtual scitbx::af::versa<double, scitbx::af::c_grid<2> >
    get_B_covariance_at_scan_point(std::size_t index) const = 0;
    virtual scitbx::af::versa<double, scitbx::af::c_grid<3> >
    get_B_covariance_at_scan_points() const = 0;
    // Get the cell parameter standard deviation
    virtual scitbx::af::small<double, 6> get_cell_parameter_sd() = 0;
    virtual scitbx::af::small<double, 6> get_cell_parameter_sd_no_calc() const = 0;
    virtual scitbx::af::small<double, 6> get_cell_parameter_sd_at_scan_point(
      std::size_t index) = 0;
    // Get the cell volume standard deviation
    virtual double get_cell_volume_sd_no_calc() const = 0;
    // Get the cell volume standard deviation
    virtual double get_cell_volume_sd() = 0;
    // Get the recalculated cell volume standard deviation
    virtual double get_recalculated_cell_volume_sd() const = 0;
    // Set the recalculated cell volume standard deviation
    virtual void set_recalculated_cell_volume_sd(
      double recalculated_cell_volume_sd) = 0;
    virtual void calc_cell_parameter_sd() = 0;
    // Reset unit cell errors
    virtual void reset_unit_cell_errors() = 0;
    // Access recalculated unit cell, intended for use by post-integration
    // refinement algorithms, such as that of dials.two_theta_refine
    virtual void set_recalculated_unit_cell(
      const cctbx::uctbx::unit_cell &unit_cell) = 0;
    virtual boost::optional<cctbx::uctbx::unit_cell> get_recalculated_unit_cell()
      const = 0;
    virtual void set_recalculated_cell_parameter_sd(
      const scitbx::af::small<double, 6> &unit_cell_sd) = 0;
    virtual scitbx::af::small<double, 6> get_recalculated_cell_parameter_sd() const = 0;
  };

  /**
   * Simple model for the crystal lattice geometry and symmetry
   *
   * A crystal is initialised from the elements of its real space axes
   * a, b, and c. Space group information must also be provided, either
   * in the form of a symbol, or an existing
   * cctbx.sgtbx.space_group object. If space_group_symbol is provided,
   * it is passed to the cctbx.sgtbx.space_group_symbols constructor.
   * This accepts either extended Hermann Mauguin format, or Hall format
   * with the prefix 'Hall:'. E.g.

   * space_group_symbol = "P b a n:1"
   *     or
   * space_group_symbol = "Hall:P 2 2 -1ab"

   */
  class Crystal : public CrystalBase {
  public:
    /**
     * Copy crystal model
     */
    Crystal(const Crystal &other)
        : space_group_(other.space_group_),
          unit_cell_(other.unit_cell_),
          recalculated_unit_cell_(other.recalculated_unit_cell_),
          U_(other.U_),
          B_(other.B_),
          A_at_scan_points_(other.A_at_scan_points_.begin(),
                            other.A_at_scan_points_.end()),
          cov_B_(other.cov_B_.accessor()),
          cell_sd_(other.cell_sd_),
          recalculated_cell_sd_(other.recalculated_cell_sd_),
          cell_volume_sd_(other.cell_volume_sd_),
          recalculated_cell_volume_sd_(other.recalculated_cell_volume_sd_) {
      std::copy(other.cov_B_.begin(), other.cov_B_.end(), cov_B_.begin());
    }

    /**
     * Initialise the crystal
     *
     * @param real_space_a The real_space a axis
     * @param real_space_b The real_space b axis
     * @param real_space_c The real_space c axis
     * @param space_group The space group object
     */
    Crystal(const vec3<double> &real_space_a,
            const vec3<double> &real_space_b,
            const vec3<double> &real_space_c,
            const cctbx::sgtbx::space_group &space_group)
        : space_group_(space_group),
          recalculated_unit_cell_(boost::none),
          cell_volume_sd_(0),
          recalculated_cell_volume_sd_(0) {
      // Setting matrix at initialisation
      mat3<double> A = mat3<double>(real_space_a[0],
                                    real_space_a[1],
                                    real_space_a[2],
                                    real_space_b[0],
                                    real_space_b[1],
                                    real_space_b[2],
                                    real_space_c[0],
                                    real_space_c[1],
                                    real_space_c[2])
                         .inverse();

      // unit cell
      set_unit_cell(real_space_a, real_space_b, real_space_c);

      // reciprocal space orthogonalisation matrix, i.e. the transpose of the
      // real space fractionalisation matrix.
      // see lecture notes from Kevin Cowtan, 'Coordinate systems, operators,
      // and transformations'
      // https://www.iucr.org/resources/commissions/computing/schools/siena-2005-crystallographic-computing-school/speakers-notes
      update_B();

      // initial orientation matrix
      U_ = A * B_.inverse();

      // set up attributes for scan-varying model
      reset_scan_points();

      // set up attributes for unit cell errors
      reset_unit_cell_errors();
    }

    /**
     * Initialise the crystal
     *
     * @param A The A matrix
     * @param space_group The space group object
     * @param reciprocal Whether or not the matrix is the direct or reciprocal
     *         matrix (default: reciprocal)
     */
    Crystal(const mat3<double> &A,
            const cctbx::sgtbx::space_group &space_group,
            const bool &reciprocal = true)
        : space_group_(space_group),
          cell_volume_sd_(0),
          recalculated_cell_volume_sd_(0) {
      if (reciprocal) {
        set_A(A);
      } else {
        set_A(A.inverse());
      }
    }

    /**
     * Constructor for pickling
     */
    Crystal(cctbx::sgtbx::space_group space_group,
            cctbx::uctbx::unit_cell unit_cell,
            boost::optional<cctbx::uctbx::unit_cell> recalculated_unit_cell,
            mat3<double> U,
            mat3<double> B,
            scitbx::af::shared<mat3<double> > A_at_scan_points,
            scitbx::af::versa<double, scitbx::af::c_grid<2> > cov_B,
            scitbx::af::small<double, 6> cell_sd,
            scitbx::af::small<double, 6> recalculated_cell_sd,
            double cell_volume_sd,
            double recalculated_cell_volume_sd)
        : space_group_(space_group),
          unit_cell_(unit_cell),
          recalculated_unit_cell_(recalculated_unit_cell),
          U_(U),
          B_(B),
          A_at_scan_points_(A_at_scan_points),
          cov_B_(cov_B),
          cell_sd_(cell_sd),
          recalculated_cell_sd_(recalculated_cell_sd),
          cell_volume_sd_(cell_volume_sd),
          recalculated_cell_volume_sd_(recalculated_cell_volume_sd) {}

    /**
     * Set the unit cell parameters
     *
     * @param real_space_a The real_space a axis
     * @param real_space_b The real_space b axis
     * @param real_space_c The real_space c axis
     */
    void set_unit_cell(const vec3<double> &real_space_a,
                       const vec3<double> &real_space_b,
                       const vec3<double> &real_space_c) {
      unit_cell_ =
        cctbx::uctbx::unit_cell(double6(real_space_a.length(),
                                        real_space_b.length(),
                                        real_space_c.length(),
                                        rad_as_deg(real_space_b.angle(real_space_c)),
                                        rad_as_deg(real_space_c.angle(real_space_a)),
                                        rad_as_deg(real_space_a.angle(real_space_b))));
      update_B();
    }

    /**
     * Set the unit cell parameters
     *
     * @param unit_cell The updated unit cell
     */
    void set_unit_cell(const cctbx::uctbx::unit_cell &unit_cell) {
      unit_cell_ = unit_cell;
      update_B();
    }

    /**
     * Update the B matrix
     */
    void update_B() {
      B_ = unit_cell_.fractionalization_matrix().transpose();
    }

    /**
     * Set the U matrix
     *
     * @param U The input matrix
     */
    void set_U(const mat3<double> &U) {
      DXTBX_ASSERT(detail::is_r3_rotation_matrix(U));
      U_ = U;
      reset_scan_points();
    }

    /**
     * Set the B matrix
     *
     * @param B The input matrix
     */
    void set_B(const mat3<double> &B) {
      // also set the unit cell
      cctbx::crystal_orientation co(B, true);
      unit_cell_ = co.unit_cell();
      B_ = unit_cell_.fractionalization_matrix().transpose();

      // reset scan-varying data, if the static B has changed
      reset_scan_points();

      // reset unit cell errors
      reset_unit_cell_errors();
    }

    /**
     * Set the A matrix
     *
     * @param A The input matrix
     */
    void set_A(const mat3<double> &A) {
      cctbx::uctbx::unit_cell uc(A.transpose().inverse());
      mat3<double> B = uc.fractionalization_matrix().transpose();
      mat3<double> U = A * B.inverse();
      set_B(B);
      set_U(U);
    }

    /**
     * @returns the U matrix
     */
    mat3<double> get_U() const {
      return U_;
    }

    /**
     * @returns the B matrix
     */
    mat3<double> get_B() const {
      return B_;
    }

    /**
     * @returns the A matrix
     */
    mat3<double> get_A() const {
      return U_ * B_;
    }

    /**
     * @returns the unit_cell
     */
    cctbx::uctbx::unit_cell get_unit_cell() const {
      return unit_cell_;
    }

    /**
     * @returns the real space vectors
     */
    scitbx::af::shared<vec3<double> > get_real_space_vectors() const {
      mat3<double> A_inv = get_A().inverse();
      scitbx::af::shared<vec3<double> > real_space_vectors;
      real_space_vectors.push_back(vec3<double>(A_inv[0], A_inv[1], A_inv[2]));
      real_space_vectors.push_back(vec3<double>(A_inv[3], A_inv[4], A_inv[5]));
      real_space_vectors.push_back(vec3<double>(A_inv[6], A_inv[7], A_inv[8]));
      return real_space_vectors;
    }

    /**
     * Set the space group
     */
    void set_space_group(cctbx::sgtbx::space_group space_group) {
      space_group_ = space_group;
    }

    /**
     * Get the space group
     */
    cctbx::sgtbx::space_group get_space_group() const {
      return space_group_;
    }

    /**
     * @returns the number of scan points
     */
    std::size_t get_num_scan_points() const {
      return A_at_scan_points_.size();
    }

    /**
     * Set the A matrix at scan points
     */
    void set_A_at_scan_points(const scitbx::af::const_ref<mat3<double> > &A) {
      A_at_scan_points_ = scitbx::af::shared<mat3<double> >(A.begin(), A.end());
    }

    /**
     * Get the A matrix at scan points
     */
    scitbx::af::shared<mat3<double> > get_A_at_scan_points() const {
      return A_at_scan_points_;
    }

    /**
     * Get the A matrix at the scan point
     */
    mat3<double> get_A_at_scan_point(std::size_t index) const {
      DXTBX_ASSERT(index < A_at_scan_points_.size());
      return A_at_scan_points_[index];
    }

    /**
     * Get the U matrix at the scan point
     */
    mat3<double> get_U_at_scan_point(std::size_t index) const {
      mat3<double> A = get_A_at_scan_point(index);
      mat3<double> B = get_B_at_scan_point(index);
      return A * B.inverse();
    }

    /**
     * Get the B matrix at the scan point
     */
    mat3<double> get_B_at_scan_point(std::size_t index) const {
      return get_unit_cell_at_scan_point(index).fractionalization_matrix().transpose();
    }

    /**
     * Get the unit cell at the scan point
     */
    cctbx::uctbx::unit_cell get_unit_cell_at_scan_point(std::size_t index) const {
      mat3<double> A = get_A_at_scan_point(index);
      return cctbx::uctbx::unit_cell(A.transpose().inverse());
    }

    /**
     * Reset the scan points
     */
    void reset_scan_points() {
      A_at_scan_points_.clear();
      cov_B_at_scan_points_ = scitbx::af::versa<double, scitbx::af::c_grid<3> >();
    }

    /**
     * Returns a copy of the current crystal model transformed by the given
     * change of basis operator to the new basis.

     * @param change_of_basis_op The change of basis operator.
     * @returns The crystal model transformed to the new basis.
     */
    boost::shared_ptr<CrystalBase> change_basis(
      cctbx::sgtbx::change_of_basis_op change_of_basis_op) const {
      // cctbx change of basis matrices and those Giacovazzo are related by
      // inverse and transpose, i.e. Giacovazzo's "M" is related to the cctbx
      // cb_op as follows:
      //   M = cb_op.c_inv().r().transpose()
      //   M_inverse = cb_op_to_minimum.c().r().transpose()

      // (Giacovazzo calls the direct matrix "A",
      //  we call the reciprocal matrix "A")
      // Therefore, from equation 2.19 in Giacovazzo:
      //   A' = M A

      // and:
      //   (A')^-1 = (M A)^-1
      //   (A')^-1 = A^-1 M^-1

      mat3<double> direct_matrix = get_A().inverse();
      mat3<double> M = change_of_basis_op.c_inv().r().transpose().as_double();

      // equation 2.19 of Giacovazzo
      mat3<double> new_direct_matrix = M * direct_matrix;
      vec3<double> real_space_a(
        new_direct_matrix[0], new_direct_matrix[1], new_direct_matrix[2]);
      vec3<double> real_space_b(
        new_direct_matrix[3], new_direct_matrix[4], new_direct_matrix[5]);
      vec3<double> real_space_c(
        new_direct_matrix[6], new_direct_matrix[7], new_direct_matrix[8]);
      // FIXME use a copy constructor. As written, this doesn't copy mosaicity or
      // parameters from derived classes
      boost::shared_ptr<Crystal> other = boost::shared_ptr<Crystal>(
        new Crystal(real_space_a,
                    real_space_b,
                    real_space_c,
                    get_space_group().change_basis(change_of_basis_op)));
      if (get_num_scan_points() > 0) {
        mat3<double> M_inv = M.inverse();
        scitbx::af::shared<mat3<double> > new_A_at_scan_points;
        for (std::size_t i = 0; i < get_num_scan_points(); ++i) {
          new_A_at_scan_points.push_back(get_A_at_scan_point(i) * M_inv);
        }
        other->set_A_at_scan_points(new_A_at_scan_points.const_ref());
      }
      if (recalculated_unit_cell_) {
        other->set_recalculated_unit_cell(
          recalculated_unit_cell_->change_basis(change_of_basis_op));
      }
      return other;
    }

    /**
     * Update crystal with parameters from another model
     */
    void update(const CrystalBase &other) {
      (*this) = *dynamic_cast<const Crystal *>(&other);
    }

    /**
     * Rotate the model around an axis and angle
     *
     * @param axis The axis to rotate around
     * @param angle The angle to rotate around
     * @param deg Degrees or radians
     */
    void rotate_around_origin(vec3<double> axis, double angle, bool deg = false) {
      // Create the rotation matrix
      mat3<double> R =
        scitbx::math::r3_rotation::axis_and_angle_as_matrix(axis, angle, deg);

      // Update U
      U_ = R * U_;

      // Update A at scan points
      for (std::size_t i = 0; i < get_num_scan_points(); ++i) {
        mat3<double> At = get_A_at_scan_point(i);
        mat3<double> Bt = get_B_at_scan_point(i);
        mat3<double> Ut = R * (At * Bt.inverse());
        A_at_scan_points_[i] = Ut * Bt;
      }
    }

    /**
     * Check if the models are equal
     */
    bool operator==(const CrystalBase &other) const {
      double d_U = 0.0;
      double d_B = 0.0;
      double eps = 1e-7;
      for (std::size_t i = 0; i < 9; ++i) {
        d_U += std::abs(U_[i] - other.get_U()[i]);
        d_B += std::abs(B_[i] - other.get_B()[i]);
      }
      if (get_num_scan_points() > 0) {
        if (get_num_scan_points() != other.get_num_scan_points()) {
          return false;
        }
        for (std::size_t j = 0; j < get_num_scan_points(); ++j) {
          mat3<double> A1 = get_A_at_scan_point(j);
          mat3<double> A2 = other.get_A_at_scan_point(j);
          double d_A = 0.0;
          for (std::size_t i = 0; i < 9; ++i) {
            d_A += std::abs(A1[i] - A2[i]);
          }
          if (d_A > eps) {
            return false;
          }
        }
      }
      return (d_U <= eps && d_B <= eps && space_group_ == other.get_space_group());
    }

    /**
     * Check if models are similar
     */
    bool is_similar_to(const CrystalBase &other,
                       double angle_tolerance = 0.01,
                       double uc_rel_length_tolerance = 0.01,
                       double uc_abs_angle_tolerance = 1.0) const {
      // space group test
      if (get_space_group() != other.get_space_group()) {
        return false;
      }

      // static orientation test
      mat3<double> U_a = get_U();
      mat3<double> U_b = other.get_U();
      DXTBX_ASSERT(detail::is_r3_rotation_matrix(U_a));
      DXTBX_ASSERT(detail::is_r3_rotation_matrix(U_b));
      mat3<double> R_ab = U_b * U_a.transpose();
      scitbx::af::tiny<double, 4> uq =
        scitbx::math::r3_rotation::matrix_as_unit_quaternion(R_ab);
      // https://github.com/cctbx/cctbx_project/issues/43 - more stable
      // mapping from unit quaternion to angle
      double x = std::sqrt(uq[1] * uq[1] + uq[2] * uq[2] + uq[3] * uq[3]);
      double angle = rad_as_deg(2.0 * std::atan2(x, uq[0]));
      if (std::abs(angle) > angle_tolerance) {
        return false;
      }

      // static unit cell test
      cctbx::uctbx::unit_cell uc_a = get_unit_cell();
      cctbx::uctbx::unit_cell uc_b = other.get_unit_cell();
      if (!uc_a.is_similar_to(uc_b, uc_rel_length_tolerance, uc_abs_angle_tolerance)) {
        return false;
      }

      // recalculated unit cell test, if both exist
      boost::optional<cctbx::uctbx::unit_cell> recalculated_uc_a =
        get_recalculated_unit_cell();
      boost::optional<cctbx::uctbx::unit_cell> recalculated_uc_b =
        other.get_recalculated_unit_cell();
      if (recalculated_uc_a && recalculated_uc_b) {
        if (!recalculated_uc_a->is_similar_to(
              *recalculated_uc_b, uc_rel_length_tolerance, uc_abs_angle_tolerance)) {
          return false;
        }
      }

      // scan varying tests
      if (get_num_scan_points() != other.get_num_scan_points()) {
        return false;
      }
      for (std::size_t i = 0; i < get_num_scan_points(); ++i) {
        mat3<double> Ut_a = get_U_at_scan_point(i);
        mat3<double> Ut_b = other.get_U_at_scan_point(i);
        DXTBX_ASSERT(detail::is_r3_rotation_matrix(Ut_a));
        DXTBX_ASSERT(detail::is_r3_rotation_matrix(Ut_b));
        mat3<double> Rt_ab = Ut_b * Ut_a.transpose();
        scitbx::af::tiny<double, 4> uq =
          scitbx::math::r3_rotation::matrix_as_unit_quaternion(Rt_ab);
        double angle = rad_as_deg(2.0 * std::acos(uq[0]));
        if (angle > angle_tolerance) {
          return false;
        }
        mat3<double> Bt_a = get_B_at_scan_point(i);
        mat3<double> Bt_b = other.get_B_at_scan_point(i);
        cctbx::uctbx::unit_cell uc_a =
          cctbx::crystal_orientation(Bt_a, true).unit_cell();
        cctbx::uctbx::unit_cell uc_b =
          cctbx::crystal_orientation(Bt_b, true).unit_cell();
        if (!uc_a.is_similar_to(
              uc_b, uc_rel_length_tolerance, uc_abs_angle_tolerance)) {
          return false;
        }
      }

      return true;
    }

    /** Check if models are not equal */
    bool operator!=(const CrystalBase &rhs) const {
      return !(*this == rhs);
    }

    /**
     * Return the 9*9 covariance matrix of elements of B, if available,
     * otherwise None. The order of elements of this matrix are determined by
     * 'flattening' B to form a 1D vector using row major ordering.
     */
    scitbx::af::versa<double, scitbx::af::c_grid<2> > get_B_covariance() const {
      return cov_B_;
    }

    /**
     * Set the covariance matrix
     */
    void set_B_covariance(
      const scitbx::af::const_ref<double, scitbx::af::c_grid<2> > &cov) {
      if (cov.accessor()[0] == 0 && cov.accessor()[1] == 0) {
        cov_B_ = scitbx::af::versa<double, scitbx::af::c_grid<2> >(cov.accessor());
      } else {
        DXTBX_ASSERT(cov.accessor()[0] == 9);
        DXTBX_ASSERT(cov.accessor()[1] == 9);
        cov_B_ = scitbx::af::versa<double, scitbx::af::c_grid<2> >(cov.accessor());
        std::copy(cov.begin(), cov.end(), cov_B_.begin());
      }
    }

    /**
     * Set the covariance matrix at scan points. The setting matrix A should be
     * set at scan-points first, so that the number of scan points is already
     * known.
     */
    void set_B_covariance_at_scan_points(
      const scitbx::af::const_ref<double, scitbx::af::c_grid<3> > &cov) {
      if (cov.accessor()[0] == 0) {
        return;
      }
      DXTBX_ASSERT(cov.accessor()[0] == get_num_scan_points());
      DXTBX_ASSERT(cov.accessor()[1] == 9);
      DXTBX_ASSERT(cov.accessor()[2] == 9);
      cov_B_at_scan_points_ =
        scitbx::af::versa<double, scitbx::af::c_grid<3> >(cov.accessor());
      std::copy(cov.begin(), cov.end(), cov_B_at_scan_points_.begin());
    }

    scitbx::af::versa<double, scitbx::af::c_grid<2> > get_B_covariance_at_scan_point(
      std::size_t index) const {
      DXTBX_ASSERT(index < cov_B_at_scan_points_.accessor()[0]);

      // Get a const ref to just the data for the scan point
      scitbx::af::const_ref<double> slice(&cov_B_at_scan_points_[index * 81], 9 * 9);

      // Copy the data and return
      scitbx::af::versa<double, scitbx::af::c_grid<2> > result(
        scitbx::af::c_grid<2>(9, 9));
      std::copy(slice.begin(), slice.end(), result.begin());
      return result;
    }

    scitbx::af::versa<double, scitbx::af::c_grid<3> > get_B_covariance_at_scan_points()
      const {
      return cov_B_at_scan_points_;
    }

    /**
     * Get the cell parameter standard deviation
     */
    scitbx::af::small<double, 6> get_cell_parameter_sd() {
      if (cov_B_.size() == 0) {
        return scitbx::af::small<double, 6>();
      } else if (cell_sd_.size() != 6) {
        calc_cell_parameter_sd();
      }
      return cell_sd_;
    }

    scitbx::af::small<double, 6> get_cell_parameter_sd_at_scan_point(
      std::size_t index) {
      mat3<double> B = get_B_at_scan_point(index);
      scitbx::af::versa<double, scitbx::af::c_grid<2> > cov_B =
        get_B_covariance_at_scan_point(index);
      scitbx::af::small<double, 6> cell_sd;
      double cell_volume;
      calc_cell_parameter_sd(B, cov_B, cell_sd, cell_volume);
      return cell_sd;
    }

    /**
     * Get the cell parameter standard deviation
     */
    scitbx::af::small<double, 6> get_cell_parameter_sd_no_calc() const {
      return cell_sd_;
    }

    /**
     * Get the cell volume standard deviation
     */
    double get_cell_volume_sd_no_calc() const {
      return cell_volume_sd_;
    }

    /**
     * Get the cell volume standard deviation
     */
    double get_recalculated_cell_volume_sd() const {
      return recalculated_cell_volume_sd_;
    }

    /**
     * Get the cell volume standard deviation
     */
    double get_cell_volume_sd() {
      if (cov_B_.size() == 0) {
        return -1;
      } else if (cell_volume_sd_ <= 0) {
        calc_cell_parameter_sd();
      }
      return cell_volume_sd_;
    }

    void calc_cell_parameter_sd() {
      calc_cell_parameter_sd(B_, cov_B_, cell_sd_, cell_volume_sd_);
    }

    void calc_cell_parameter_sd(
      const mat3<double> &B,
      const scitbx::af::versa<double, scitbx::af::c_grid<2> > &cov_B,
      scitbx::af::small<double, 6> &cell_sd,
      double &cell_volume_sd) {
      // self._cov_B is the covariance matrix of elements of the B matrix. We
      // need to construct the covariance matrix of elements of the
      // transpose of B. The vector of elements of B is related to the
      // vector of elements of its transpose by a permutation, P.
      double P[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

      // We can use P to replace var_cov with the covariance matrix of
      // elements of the transpose of B.
      scitbx::af::small<double, 81> PcovB(81, 0);
      scitbx::af::small<double, 81> var_cov(81, 0);
      scitbx::matrix::multiply(&P[0], &cov_B[0], 9, 9, 9, &PcovB[0]);
      scitbx::matrix::multiply(&PcovB[0], &P[0], 9, 9, 9, &var_cov[0]);

      // From B = (O^-1)^T we can convert this
      // to the covariance matrix of the real space orthogonalisation matrix
      mat3<double> Bt = B.transpose();
      mat3<double> O = Bt.inverse();
      scitbx::af::versa<double, scitbx::af::c_grid<2> > cov_O =
        detail::matrix_inverse_error_propagation(
          scitbx::af::ref<double, scitbx::af::c_grid<2> >(&Bt[0],
                                                          scitbx::af::c_grid<2>(3, 3)),
          scitbx::af::ref<double, scitbx::af::c_grid<2> >(&var_cov[0],
                                                          scitbx::af::c_grid<2>(9, 9)));

      // The real space unit cell vectors are given by
      vec3<double> vec_a = O * vec3<double>(1, 0, 0);
      vec3<double> vec_b = O * vec3<double>(0, 1, 0);
      vec3<double> vec_c = O * vec3<double>(0, 0, 1);

      // So the unit cell parameters are
      double a = vec_a.length();
      double b = vec_b.length();
      double c = vec_c.length();
      /* double alpha = std::acos((vec_b * vec_c) / (b*c)); */
      /* double beta  = std::acos((vec_a * vec_c) / (a*c)); */
      /* double gamma = std::acos((vec_a * vec_b) / (a*b)); */

      /*
       * The estimated errors are calculated by error propagation from cov_O. In
       * each case we define a function F(O) that converts the matrix O into the
       * unit cell parameter of interest. To do error propagation to get the
       * variance of that cell parameter we need the Jacobian of the function.
       * This is a 1*9 matrix of derivatives in the order of the elements of O
       *
       *     / dF   dF   dF   dF   dF   dF   dF   dF   dF  \
       * J = | ---, ---, ---, ---, ---, ---, ---, ---, --- |
       *     \ da1  db1  dc1  da2  db2  dc2  da3  db3  dc3 /
       *
       */

      // For cell length |a|, F = sqrt(a1^2 + a2^2 + a3^2)
      double jacobian1[9] = {
        vec_a[0] / a, 0, 0, vec_a[1] / a, 0, 0, vec_a[2] / a, 0, 0};
      double var_a = detail::multiply_A_B_AT<double>(jacobian1, &cov_O[0], 1, 9, 9);

      // For cell length |b|, F = sqrt(b1^2 + b2^2 + b3^2)
      double jacobian2[9] = {
        0, vec_b[0] / b, 0, 0, vec_b[1] / b, 0, 0, vec_b[2] / b, 0};
      double var_b = detail::multiply_A_B_AT<double>(jacobian2, &cov_O[0], 1, 9, 9);

      // For cell length |c|, F = sqrt(c1^2 + c2^2 + c3^2)
      double jacobian3[9] = {
        0, 0, vec_c[0] / c, 0, 0, vec_c[1] / c, 0, 0, vec_c[2] / c};
      double var_c = detail::multiply_A_B_AT<double>(jacobian3, &cov_O[0], 1, 9, 9);

      // For cell volume (a X b).c,
      // F = c1(a2*b3 - b2*a3) + c2(a3*b1 - b3*a1) + c3(a1*b2 - b1*a2)
      double a1 = vec_a[0];
      double a2 = vec_a[1];
      double a3 = vec_a[2];
      double b1 = vec_b[0];
      double b2 = vec_b[1];
      double b3 = vec_b[2];
      double c1 = vec_c[0];
      double c2 = vec_c[1];
      double c3 = vec_c[2];
      double jacobian4[9] = {c3 * b2 - c2 * b3,
                             c1 * b3 - c3 * b1,
                             c2 * b1 - c1 * b2,
                             c2 * a3 - c3 * a2,
                             c3 * a1 - c1 * a3,
                             c1 * a2 - c2 * a1,
                             a2 * b3 - b2 * a3,
                             a3 * b1 - b3 * a1,
                             a1 * b2 - b1 * a2};
      double var_V = detail::multiply_A_B_AT<double>(jacobian4, &cov_O[0], 1, 9, 9);
      cell_volume_sd = std::sqrt(var_V);

      // For the unit cell angles we need to calculate derivatives of the angles
      // with respect to the elements of O
      scitbx::af::tiny<vec3<double>, 2> dalpha_dbdc =
        scitbx::math::angle_derivative_wrt_vectors(vec_b, vec_c);
      scitbx::af::tiny<vec3<double>, 2> dbeta_dadc =
        scitbx::math::angle_derivative_wrt_vectors(vec_a, vec_c);
      scitbx::af::tiny<vec3<double>, 2> dgamma_dadb =
        scitbx::math::angle_derivative_wrt_vectors(vec_a, vec_b);
      vec3<double> dalpha_db = dalpha_dbdc[0];
      vec3<double> dalpha_dc = dalpha_dbdc[1];
      vec3<double> dbeta_da = dbeta_dadc[0];
      vec3<double> dbeta_dc = dbeta_dadc[1];
      vec3<double> dgamma_da = dgamma_dadb[0];
      vec3<double> dgamma_db = dgamma_dadb[1];

      // For angle alpha, F = acos( b.c / |b||c|)
      double jacobian5[9] = {0,
                             dalpha_db[0],
                             dalpha_dc[0],
                             0,
                             dalpha_db[1],
                             dalpha_dc[1],
                             0,
                             dalpha_db[2],
                             dalpha_dc[2]};
      double var_alpha = detail::multiply_A_B_AT<double>(jacobian5, &cov_O[0], 1, 9, 9);

      // For angle beta, F = acos( a.c / |a||c|)
      double jacobian6[9] = {dbeta_da[0],
                             0,
                             dbeta_dc[0],
                             dbeta_da[1],
                             0,
                             dbeta_dc[1],
                             dbeta_da[2],
                             0,
                             dbeta_dc[2]};
      double var_beta = detail::multiply_A_B_AT<double>(jacobian6, &cov_O[0], 1, 9, 9);

      // For angle gamma, F = acos( a.b / |a||b|)
      double jacobian7[9] = {dgamma_da[0],
                             dgamma_db[0],
                             0,
                             dgamma_da[1],
                             dgamma_db[1],
                             0,
                             dgamma_da[2],
                             dgamma_db[2],
                             0};
      double var_gamma = detail::multiply_A_B_AT<double>(jacobian7, &cov_O[0], 1, 9, 9);

      // Symmetry constraints may mean variances of the angles should be zero.
      // Floating point error could lead to negative variances. Ensure these are
      // caught before taking their square root!
      var_alpha = std::max(0.0, var_alpha);
      var_beta = std::max(0.0, var_beta);
      var_gamma = std::max(0.0, var_gamma);

      // Set the cell sd values
      cell_sd.resize(6);
      cell_sd[0] = std::sqrt(var_a);
      cell_sd[1] = std::sqrt(var_b);
      cell_sd[2] = std::sqrt(var_c);
      cell_sd[3] = rad_as_deg(std::sqrt(var_alpha));
      cell_sd[4] = rad_as_deg(std::sqrt(var_beta));
      cell_sd[5] = rad_as_deg(std::sqrt(var_gamma));
    }

    /**
     * Reset unit cell errors
     */
    void reset_unit_cell_errors() {
      cov_B_ = scitbx::af::versa<double, scitbx::af::c_grid<2> >();
      cov_B_at_scan_points_ = scitbx::af::versa<double, scitbx::af::c_grid<3> >();
      cell_sd_ = scitbx::af::small<double, 6>();
      cell_volume_sd_ = 0;
      recalculated_cell_sd_ = scitbx::af::small<double, 6>();
      recalculated_cell_volume_sd_ = 0;
    }

    void set_recalculated_unit_cell(const cctbx::uctbx::unit_cell &unit_cell) {
      recalculated_unit_cell_ = unit_cell;
      recalculated_cell_sd_ = scitbx::af::small<double, 6>();
      recalculated_cell_volume_sd_ = 0;
    }

    boost::optional<cctbx::uctbx::unit_cell> get_recalculated_unit_cell() const {
      return recalculated_unit_cell_;
    }

    void set_recalculated_cell_parameter_sd(
      const scitbx::af::small<double, 6> &unit_cell_sd) {
      recalculated_cell_sd_ = unit_cell_sd;
    }

    void set_recalculated_cell_volume_sd(double recalculated_cell_volume_sd) {
      recalculated_cell_volume_sd_ = recalculated_cell_volume_sd;
    }

    scitbx::af::small<double, 6> get_recalculated_cell_parameter_sd() const {
      return recalculated_cell_sd_;
    }

  protected:
    cctbx::sgtbx::space_group space_group_;
    cctbx::uctbx::unit_cell unit_cell_;
    boost::optional<cctbx::uctbx::unit_cell> recalculated_unit_cell_;
    mat3<double> U_;
    mat3<double> B_;
    scitbx::af::shared<mat3<double> > A_at_scan_points_;
    scitbx::af::versa<double, scitbx::af::c_grid<2> > cov_B_;
    scitbx::af::versa<double, scitbx::af::c_grid<3> > cov_B_at_scan_points_;
    scitbx::af::small<double, 6> cell_sd_;
    scitbx::af::small<double, 6> recalculated_cell_sd_;
    double cell_volume_sd_;
    double recalculated_cell_volume_sd_;
  };

  /* Extended Crystal class adding a simple value for mosaicity.
   * The deg parameter controls whether this value is treated as being an
   * angle in degrees or radians.
   */
  class MosaicCrystalKabsch2010 : public Crystal {
  public:
    /**
     * Copy crystal model from MosaicCrystalKabsch2010 class
     */
    MosaicCrystalKabsch2010(const MosaicCrystalKabsch2010 &other)
        : Crystal(other), mosaicity_(other.mosaicity_) {}
    /**
     * Copy crystal model from Crystal class
     */
    MosaicCrystalKabsch2010(const Crystal &other) : Crystal(other), mosaicity_(0) {}

    /**
     * Initialise the crystal
     *
     * @param real_space_a The real_space a axis
     * @param real_space_b The real_space b axis
     * @param real_space_c The real_space c axis
     * @param space_group The space group object
     */
    MosaicCrystalKabsch2010(const vec3<double> &real_space_a,
                            const vec3<double> &real_space_b,
                            const vec3<double> &real_space_c,
                            const cctbx::sgtbx::space_group &space_group)
        : Crystal(real_space_a, real_space_b, real_space_c, space_group),
          mosaicity_(0) {}

    /**
     * Constructor for pickling
     */
    MosaicCrystalKabsch2010(
      cctbx::sgtbx::space_group space_group,
      cctbx::uctbx::unit_cell unit_cell,
      boost::optional<cctbx::uctbx::unit_cell> recalculated_unit_cell,
      mat3<double> U,
      mat3<double> B,
      scitbx::af::shared<mat3<double> > A_at_scan_points,
      scitbx::af::versa<double, scitbx::af::c_grid<2> > cov_B,
      scitbx::af::small<double, 6> cell_sd,
      scitbx::af::small<double, 6> recalculated_cell_sd,
      double cell_volume_sd,
      double recalculated_cell_volume_sd,
      double mosaicity)
        : Crystal(space_group,
                  unit_cell,
                  recalculated_unit_cell,
                  U,
                  B,
                  A_at_scan_points,
                  cov_B,
                  cell_sd,
                  recalculated_cell_sd,
                  cell_volume_sd,
                  recalculated_cell_volume_sd),
          mosaicity_(mosaicity) {}

    /**
     * Check if the models are equal
     */
    bool operator==(const CrystalBase &other) const {
      double eps = 1e-7;
      const MosaicCrystalKabsch2010 *mosaic_other =
        dynamic_cast<const MosaicCrystalKabsch2010 *>(&other);
      if (!mosaic_other) return false;

      double d_mosaicity = std::abs(mosaicity_ - mosaic_other->get_mosaicity());
      return (d_mosaicity <= eps && Crystal::operator==(other));
    }

    /**
     * Check if models are similar
     */
    bool is_similar_to(const CrystalBase &other,
                       double angle_tolerance,
                       double uc_rel_length_tolerance,
                       double uc_abs_angle_tolerance,
                       double mosaicity_tolerance) const {
      // mosaicity test
      const MosaicCrystalKabsch2010 *mosaic_other =
        dynamic_cast<const MosaicCrystalKabsch2010 *>(&other);
      if (!mosaic_other) return false;

      double m_a = get_mosaicity();
      double m_b = mosaic_other->get_mosaicity();
      double min_m = std::min(m_a, m_b);
      double max_m = std::max(m_a, m_b);
      if (min_m <= 0) {
        if (max_m > 0.0) {
          return false;
        }
      } else if (min_m / max_m < mosaicity_tolerance) {
        return false;
      }

      return Crystal::is_similar_to(
        other, angle_tolerance, uc_rel_length_tolerance, uc_abs_angle_tolerance);
    }

    bool is_similar_to(const CrystalBase &other,
                       double angle_tolerance = 0.01,
                       double uc_rel_length_tolerance = 0.01,
                       double uc_abs_angle_tolerance = 1.0) const /* override */ {
      // Call with a default mosaicity tolerance
      return is_similar_to(
        other, angle_tolerance, uc_rel_length_tolerance, uc_abs_angle_tolerance, 0.8);
    }

    /**
     * Get the mosaicity
     */
    double get_mosaicity(bool deg = true) const {
      if (deg == false) {
        return deg_as_rad(mosaicity_);
      }
      return mosaicity_;
    }

    /**
     * Set the mosaicity
     */
    void set_mosaicity(double mosaicity, bool deg = true) {
      if (deg == true) {
        mosaicity_ = mosaicity;
      } else {
        mosaicity_ = rad_as_deg(mosaicity_);
      }
    }

    double mosaicity_;
  };

  /* Extended Crystal class adding two values for mosaicity, the mosaic
   * angle and the domain size.
   */
  class MosaicCrystalSauter2014 : public Crystal {
  public:
    /**
     * Copy crystal model from MosaicCrystalSauter2014 class
     */
    MosaicCrystalSauter2014(const MosaicCrystalSauter2014 &other)
        : Crystal(other),
          half_mosaicity_deg_(other.half_mosaicity_deg_),
          domain_size_ang_(other.domain_size_ang_) {}
    /**
     * Copy crystal model from Crystal class
     */
    MosaicCrystalSauter2014(const Crystal &other)
        : Crystal(other), half_mosaicity_deg_(0), domain_size_ang_(0) {}

    /**
     * Initialise the crystal
     *
     * @param real_space_a The real_space a axis
     * @param real_space_b The real_space b axis
     * @param real_space_c The real_space c axis
     * @param space_group The space group object
     */
    MosaicCrystalSauter2014(const vec3<double> &real_space_a,
                            const vec3<double> &real_space_b,
                            const vec3<double> &real_space_c,
                            const cctbx::sgtbx::space_group &space_group)
        : Crystal(real_space_a, real_space_b, real_space_c, space_group),
          half_mosaicity_deg_(0),
          domain_size_ang_(0) {}

    /**
     * Constructor for pickling
     */
    MosaicCrystalSauter2014(
      cctbx::sgtbx::space_group space_group,
      cctbx::uctbx::unit_cell unit_cell,
      boost::optional<cctbx::uctbx::unit_cell> recalculated_unit_cell,
      mat3<double> U,
      mat3<double> B,
      scitbx::af::shared<mat3<double> > A_at_scan_points,
      scitbx::af::versa<double, scitbx::af::c_grid<2> > cov_B,
      scitbx::af::small<double, 6> cell_sd,
      scitbx::af::small<double, 6> recalculated_cell_sd,
      double cell_volume_sd,
      double recalculated_cell_volume_sd,
      double half_mosaicity_deg,
      double domain_size_ang)
        : Crystal(space_group,
                  unit_cell,
                  recalculated_unit_cell,
                  U,
                  B,
                  A_at_scan_points,
                  cov_B,
                  cell_sd,
                  recalculated_cell_sd,
                  cell_volume_sd,
                  recalculated_cell_volume_sd),
          half_mosaicity_deg_(half_mosaicity_deg),
          domain_size_ang_(domain_size_ang) {}

    /**
     * Check if the models are equal
     */
    bool operator==(const CrystalBase &other) const {
      double eps = 1e-7;
      const MosaicCrystalSauter2014 *mosaic_other =
        dynamic_cast<const MosaicCrystalSauter2014 *>(&other);
      if (!mosaic_other) return false;

      double d_half_mosaicity_deg =
        std::abs(half_mosaicity_deg_ - mosaic_other->get_half_mosaicity_deg());
      double d_domain_size_ang =
        std::abs(domain_size_ang_ - mosaic_other->get_domain_size_ang());
      return (d_half_mosaicity_deg <= eps && d_domain_size_ang <= eps
              && Crystal::operator==(other));
    }

    /**
     * Check if models are similar
     */
    bool is_similar_to(const CrystalBase &other,
                       double angle_tolerance,
                       double uc_rel_length_tolerance,
                       double uc_abs_angle_tolerance,
                       double half_mosaicity_tolerance,
                       double domain_size_tolerance = 1.0) const {
      // mosaicity test
      const MosaicCrystalSauter2014 *mosaic_other =
        dynamic_cast<const MosaicCrystalSauter2014 *>(&other);
      if (!mosaic_other) return false;

      double m_a = get_half_mosaicity_deg();
      double m_b = mosaic_other->get_half_mosaicity_deg();
      double min_m = std::min(m_a, m_b);
      double max_m = std::max(m_a, m_b);
      if (min_m <= 0) {
        if (max_m > 0.0) {
          return false;
        }
      } else if (min_m / max_m < half_mosaicity_tolerance) {
        return false;
      }

      if (std::abs(get_domain_size_ang() - mosaic_other->get_domain_size_ang())
          > domain_size_tolerance) {
        return false;
      }

      return Crystal::is_similar_to(
        other, angle_tolerance, uc_rel_length_tolerance, uc_abs_angle_tolerance);
    }

    /// Check if models are similar
    bool is_similar_to(const CrystalBase &other,
                       double angle_tolerance = 0.01,
                       double uc_rel_length_tolerance = 0.01,
                       double uc_abs_angle_tolerance = 1.0) const /* override */ {
      // Call with a default half-mosaicity tolerance
      // Use the function-default domain_size_tolerance
      return is_similar_to(
        other, angle_tolerance, uc_rel_length_tolerance, uc_abs_angle_tolerance, 0.4);
    }

    /**
     * Get the half mosaic angle value (degrees)
     */
    double get_half_mosaicity_deg() const {
      return half_mosaicity_deg_;
    }

    /**
     * Set the half mosaic angle value (degrees)
     */
    void set_half_mosaicity_deg(double half_mosaicity_deg) {
      half_mosaicity_deg_ = half_mosaicity_deg;
    }

    /**
     * Get the half mosaic angle value (degrees)
     */
    double get_domain_size_ang() const {
      return domain_size_ang_;
    }

    /**
     * Set the mosaic domain size (angstroms)
     */
    void set_domain_size_ang(double domain_size_ang) {
      domain_size_ang_ = domain_size_ang;
    }

    double half_mosaicity_deg_;
    double domain_size_ang_;
  };

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_CRYSTAL_H
