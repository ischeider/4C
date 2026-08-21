// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_CALC_EAS_HELPERS_HPP
#define FOUR_C_SOLID_ELE_CALC_EAS_HELPERS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_solid_ele_calc_lib.hpp"

#include <cstring>

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  /*!
   * @brief Solve for the inverse of a matrix and ignore any errors
   *
   * @tparam dim : matrix dimensions
   * @param matrix(in/out) : matrix to be inverted
   */
  template <unsigned int dim>
  void solve_for_inverse_ignoring_errors(Core::LinAlg::Matrix<dim, dim>& matrix)
  {
    Core::LinAlg::FixedSizeSerialDenseSolver<dim, dim, 1> solve_for_inverse;
    solve_for_inverse.set_matrix(matrix);

    solve_for_inverse.invert();
  }

  enum class EasType
  {
    soh8_easnone,
    eastype_h8_9,
    eastype_h8_21,
    eastype_sh8_7,
    eastype_sw6_1,
    eastype_q4_4,
    eastype_undefined
  };

  template <Discret::Elements::EasType eastype>
  struct EasTypeToNumEas
  {
  };
  template <>
  struct EasTypeToNumEas<Discret::Elements::EasType::eastype_h8_9>
  {
    static constexpr int num_eas = 9;
  };
  template <>
  struct EasTypeToNumEas<Discret::Elements::EasType::eastype_h8_21>
  {
    static constexpr int num_eas = 21;
  };
  template <>
  struct EasTypeToNumEas<Discret::Elements::EasType::eastype_sh8_7>
  {
    static constexpr int num_eas = 7;
  };
  template <>
  struct EasTypeToNumEas<Discret::Elements::EasType::eastype_sw6_1>
  {
    static constexpr int num_eas = 1;
  };
  template <>
  struct EasTypeToNumEas<Discret::Elements::EasType::eastype_q4_4>
  {
    static constexpr int num_eas = 4;
  };
  template <>
  struct EasTypeToNumEas<Discret::Elements::EasType::eastype_undefined>
  {
  };

  template <Core::FE::CellType celltype>
  inline static constexpr int num_nodes = Core::FE::num_nodes(celltype);

  template <Core::FE::CellType celltype>
  inline static constexpr int num_dim = Core::FE::dim<celltype>;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_str = num_dim<celltype> * (num_dim<celltype> + 1) / 2;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_dof_per_ele = num_nodes<celltype> * num_dim<celltype>;


  /// struct for EAS matrices and vectors to be stored between iterations
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  struct EasIterationData
  {
    constexpr static int num_eas = Discret::Elements::EasTypeToNumEas<eastype>::num_eas;

    /// inverse EAS matrix K_{alpha alpha}
    Core::LinAlg::Matrix<num_eas, num_eas> invKaa{Core::LinAlg::Initialization::zero};

    /// EAS matrix K_{d alpha}
    Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, num_eas> Kda{
        Core::LinAlg::Initialization::zero};

    /// EAS enhancement vector s
    Core::LinAlg::Matrix<num_eas, 1> s{Core::LinAlg::Initialization::zero};

    /// discrete enhanced strain scalars increment
    Core::LinAlg::Matrix<num_eas, 1> alpha_inc{Core::LinAlg::Initialization::zero};

    /// discrete enhanced strain scalars alpha
    Core::LinAlg::Matrix<num_eas, 1> alpha{Core::LinAlg::Initialization::zero};
  };

  template <Core::FE::CellType celltype>
  struct CentroidTransformation
  {
    // transformation matrix T0^{-T}, which maps the matrix M from parameter space to the material
    // configuration see Andelfinger et al., EAS-elements, 1993, doi: 10.1002/nme.1620360805
    Core::LinAlg::Matrix<num_str<celltype>, num_str<celltype>> T0invT{};

    // Jacobi determinant evaluated at the element centroid
    double detJ0{};
  };



  /*!
   * @brief Evaluates and returns the centroid transformation quantities, i.e., the jacobi
   * determinant at the element centroid and the transformation matrix T0^{-T}
   *
   * @tparam celltype : Cell type
   * @param nodal_coordinates(in) : reference and current coordinates of the nodes of the element
   * @return CentroidTransformation<celltype> : Jacobi determinant at the element centroid and
   * transformation matrix T0^{-T}
   */
  template <Core::FE::CellType celltype>
  Discret::Elements::CentroidTransformation<celltype> evaluate_centroid_transformation(
      const Discret::Elements::ElementNodes<celltype>& nodal_coordinates)
  {
    Discret::Elements::CentroidTransformation<celltype> centroid_transformation;

    // 1) compute jacobian at element centroid
    const Discret::Elements::JacobianMapping<celltype> jacobian_mapping_centroid =
        Discret::Elements::evaluate_jacobian_mapping_centroid(nodal_coordinates);

    centroid_transformation.detJ0 = jacobian_mapping_centroid.determinant_;

    // 2) compute matrix T0^{-T}: T0^{-T} maps the matrix M from local to global coordinates
    centroid_transformation.T0invT =
        evaluate_voigt_transformation_matrix(jacobian_mapping_centroid);

    return centroid_transformation;
  }

  /*!
   * @brief Extracts and returns the residual displacement
   *
   * @tparam celltype : Cell type
   * @param discretization(in) : reference to the discretization
   * @param lm(in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   * @return double : residual displacement or displacement increment
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1> get_displacement_increment(
      const Core::FE::Discretization& discretization, const std::vector<int>& lm)
  {
    auto residual_from_dis = discretization.get_state("residual displacement");
    const std::array residual =
        Core::FE::extract_values_as_array<Discret::Elements::num_dof_per_ele<celltype>>(
            *residual_from_dis, lm);
    Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1> displ_inc(
        Core::LinAlg::Initialization::uninitialized);
    std::memcpy(displ_inc.data(), residual.data(), sizeof(double) * residual.size());

    return displ_inc;
  }

  /*!
   * @brief Updates the enhanced strains scalar increment
   *
   * @tparam celltype, eastype
   * @param displ_inc(in) : displacement increment delta_D_{i+1}
   * @param eas_iteration_data(in) : EAS matrices and vectors from iteration i
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  void update_alpha_increment(
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1>& displ_inc,
      Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data)
  {
    // the enhanced strains scalar increment is computed to:
    // delta_alpha_{i+1} = - invKaa_{i} (s_{i} + Kad_{i} delta_D_{i+1})

    // init as enhancement vector s_{i}
    Core::LinAlg::Matrix<Discret::Elements::EasTypeToNumEas<eastype>::num_eas, 1> tmp(
        eas_iteration_data.s);

    // addition of Kad_{i} delta_D_{i+1}
    tmp.multiply_tn(1.0, eas_iteration_data.Kda, displ_inc, 1.0);

    // multiplication with (- invKaa_{i})
    eas_iteration_data.alpha_inc.multiply(-1.0, eas_iteration_data.invKaa, tmp);
  }

  /*!
   * @brief Updates the enhanced strain scalars, alpha_inc and alpha in the iteration data
   * accordingly
   *
   * @tparam celltype, eastype
   * @param eas_iteration_data(in/out) : EAS matrices and vectors
   * @param displ_increment (in) : displacement increment
   * @param step_length(in) : Length of the step (relevant for line search solvers)
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  void update_alpha(Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data,
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1>& displ_inc,
      const double step_length = 1.0)
  {
    // compute the enhanced strain scalar increment delta_alpha
    update_alpha_increment<celltype, eastype>(displ_inc, eas_iteration_data);

    // update alpha_i with the increment delta_alpha such that alpha_{i+1} = alpha_{i} +
    // delta_alpha
    eas_iteration_data.alpha.update(step_length, eas_iteration_data.alpha_inc, 1.0);
  }

  /*!
   * @brief Correct alpha in a line search step by adapting the step length
   *
   * @tparam celltype
   * @tparam eastype
   * @param eas_iteration_data (in/out) : EAS iteration data
   * @param new_step_length (in) : new step length
   * @param old_step_length (in) : old step length
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  void correct_alpha(Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data,
      const double new_step_length, const double old_step_length)
  {
    eas_iteration_data.alpha.update(
        new_step_length - old_step_length, eas_iteration_data.alpha_inc, 1.0);
  }

  /*!
   * @brief Compute the matrix M which is the element-wise matrix of the shape functions for the
   * enhanced strains in the parameter space
   *
   * @tparam celltype, eastype
   * @param xi(in) : coordinate in the parameter space
   * @return Core::LinAlg::Matrix<num_str, num_eas> : enhanced strains shape function matrix in
   * parameter space
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
      Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
  evaluate_eas_shape_functions_parameter_space(
      const Core::LinAlg::Tensor<double, Discret::Elements::num_dim<celltype>>& xi)
  {
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
        M(Core::LinAlg::Initialization::zero);

    switch (eastype)
    {
      /* easmild is the EAS interpolation of 9 modes, based on
      **            r 0 0   0 0 0 0 0 0
      **            0 s 0   0 0 0 0 0 0
      **    M =     0 0 t   0 0 0 0 0 0
      **            0 0 0   r s 0 0 0 0
      **            0 0 0   0 0 s t 0 0
      **            0 0 0   0 0 0 0 r t
      */
      case Discret::Elements::EasType::eastype_h8_9:
      {
        M(0, 0) = xi(0);
        M(1, 1) = xi(1);
        M(2, 2) = xi(2);
        M(3, 3) = xi(0);
        M(3, 4) = xi(1);
        M(4, 5) = xi(1);
        M(4, 6) = xi(2);
        M(5, 7) = xi(0);
        M(5, 8) = xi(2);

        break;
      }
      /* easfull is the EAS interpolation of 21 modes, based on
      **            r 0 0   0 0 0 0 0 0   0  0  0  0  0  0   rs rt 0  0  0  0
      **            0 s 0   0 0 0 0 0 0   0  0  0  0  0  0   0  0  rs st 0  0
      **    M =     0 0 t   0 0 0 0 0 0   0  0  0  0  0  0   0  0  0  0  rt st
      **            0 0 0   r s 0 0 0 0   rt st 0  0  0  0   0  0  0  0  0  0
      **            0 0 0   0 0 s t 0 0   0  0  rs rt 0  0   0  0  0  0  0  0
      **            0 0 0   0 0 0 0 r t   0  0  0  0  rs st  0  0  0  0  0  0
      */
      case Discret::Elements::EasType::eastype_h8_21:
      {
        M(0, 0) = xi(0);
        M(0, 15) = xi(0) * xi(1);
        M(0, 16) = xi(0) * xi(2);
        M(1, 1) = xi(1);
        M(1, 17) = xi(0) * xi(1);
        M(1, 18) = xi(1) * xi(2);
        M(2, 2) = xi(2);
        M(2, 19) = xi(0) * xi(2);
        M(2, 20) = xi(1) * xi(2);
        M(3, 3) = xi(0);
        M(3, 4) = xi(1);
        M(3, 9) = xi(0) * xi(2);
        M(3, 10) = xi(1) * xi(2);
        M(4, 5) = xi(1);
        M(4, 6) = xi(2);
        M(4, 11) = xi(0) * xi(1);
        M(4, 12) = xi(0) * xi(2);
        M(5, 7) = xi(0);
        M(5, 8) = xi(2);
        M(5, 13) = xi(0) * xi(1);
        M(5, 14) = xi(1) * xi(2);

        break;
      }
      /* eassosh8 is the EAS interpolation for the Solid-Shell with t=thickness dir.
      ** consisting of 7 modes, based on
      **            r 0 0   0 0 0  0
      **            0 s 0   0 0 0  0
      **    M =     0 0 t   0 0 rt st
      **            0 0 0   r s 0  0
      **            0 0 0   0 0 0  0
      **            0 0 0   0 0 0  0
      */
      case Discret::Elements::EasType::eastype_sh8_7:
      {
        /* eassosh8 is the EAS interpolation for the Solid-Shell with t=thickness dir.
        ** consisting of 7 modes, based on
        **            r 0 0   0 0 0  0
        **            0 s 0   0 0 0  0
        **    M =     0 0 t   0 0 rt st
        **            0 0 0   r s 0  0
        **            0 0 0   0 0 0  0
        **            0 0 0   0 0 0  0
        */
        M(0, 0) = xi(0);
        M(1, 1) = xi(1);
        M(2, 2) = xi(2);
        M(2, 5) = xi(0) * xi(2);
        M(2, 6) = xi(1) * xi(2);
        M(3, 3) = xi(0);
        M(3, 4) = xi(1);

        break;
      }
      case Discret::Elements::EasType::eastype_sw6_1:
      {
        /* eassosw6 is the EAS interpolation for the Solid-Shell with t=thickness dir.
        ** consisting of 1 modes, based on
        **            0
        **            0
        **    M =     t
        **            0
        **            0
        **            0
        */
        M(2, 0) = xi(2);

        break;
      }
      case Discret::Elements::EasType::eastype_q4_4:
      {
        /* eastype_q4_4 is the EAS interpolation for the 2d Solid quad4 element.
        ** consisting of 4 modes, based on
        **           r 0 0 0
        **    M =    0 s 0 0
        **           0 0 r s
        */
        M(0, 0) = xi(0);
        M(1, 1) = xi(1);
        M(2, 2) = xi(0);
        M(2, 3) = xi(1);

        break;
      }
      default:
        FOUR_C_THROW("unknown EAS type");
        break;
    }
    return M;
  }

  /*!
   * @brief Map the matrix M in the parameter space to Mtilde in the material configuration and
   * return Mtilde
   *
   * @tparam celltype, eastype
   * @param detJ(in) : Jacobi determinant at Gauss point
   * @param centroid_transformation(in) : transformation matrix T0^{-T} and Jacobi determinant
   * at element centroid
   * @param M(in) : matrix M in the parameter space
   * @return Core::LinAlg::Matrix<num_str, num_eas> : matrix Mtilde in the material
   * configuration
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
      Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
  map_eas_shape_functions_to_material_config(const double detJ,
      const Discret::Elements::CentroidTransformation<celltype>& centroid_transformation,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Discret::Elements::EasTypeToNumEas<eastype>::num_eas>& M)
  {
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
        Mtilde;

    // Mtilde = detJ0/detJ T0^{-T} M
    Mtilde.multiply(centroid_transformation.detJ0 / detJ, centroid_transformation.T0invT, M);

    return Mtilde;
  }

  /*!
   * @brief Evaluate the element-wise matrix of the shape functions for the enhanced strains in
   * the parameter space Mtilde. Therefore set up M (in the material configuration) and map M to
   * Mtilde via T0^{-T}.
   *
   * @tparam celltype, eastype
   * @param detJ(in) : Jacobi determinant at Gauss point
   * @param centroid_transformation(in) : transformation matrix T0^{-T} and Jacobi determinant
   * at element centroid
   * @param xi(in) : coordinate in the parameter space
   * @return Core::LinAlg::Matrix<num_str, num_eas> : matrix Mtilde in the material
   * configuration
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
      Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
  evaluate_eas_shape_functions_material_config(const double detJ,
      const Discret::Elements::CentroidTransformation<celltype>& centroid_transformation,
      const Core::LinAlg::Tensor<double, Discret::Elements::num_dim<celltype>>& xi)
  {
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
        M(evaluate_eas_shape_functions_parameter_space<celltype, eastype>(xi));
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
        Mtilde = map_eas_shape_functions_to_material_config<celltype, eastype>(
            detJ, centroid_transformation, M);
    return Mtilde;
  }

  /*!
   * @brief Evaluates the enhanced assumed Green-Lagrange strains E^{enh} = Mtilde alpha
   *
   * Background: Choose deformation gradient F as sum of displacement-based F^{u} and enhanced
   * gradient F^{enh}. Considering F_0 the deformation gradient evaluated at the element
   * centroid, F^{enh} is computed to F^{enh} = F_0^{u} Mtilde alpha.
   *
   * @tparam celltype, eastype
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param alpha(in) : enhanced strain scalars
   * @return Core::LinAlg::Matrix<num_str, 1>  : enhanced Green-Lagrange strains E^{enh}
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
  evaluate_enhanced_assumed_gl_strains(
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Discret::Elements::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const Core::LinAlg::Matrix<Discret::Elements::EasTypeToNumEas<eastype>::num_eas, 1>& alpha)
  {
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>, 1> enhanced_gl_strain(
        Core::LinAlg::Initialization::uninitialized);
    enhanced_gl_strain.multiply(Mtilde, alpha);

    if constexpr (Core::FE::dim<celltype> == 2)
    {
      // convert to stress like voigt notation
      enhanced_gl_strain(2) *= 0.5;
    }
    else
    {
      Core::LinAlg::Voigt::Strains::to_stress_like(enhanced_gl_strain, enhanced_gl_strain);
    }
    Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        gl_strain_eas_tensor;
    std::ranges::copy_n(enhanced_gl_strain.data(), Discret::Elements::num_str<celltype>,
        gl_strain_eas_tensor.data());
    return gl_strain_eas_tensor;
  }

  /*!
   * @brief Evaluate the enhanced assumed Green-Lagrange strains E^{enh}

   * @tparam celltype, eastype
   * @param displacement_based_mapping(in) : displacement-based spatial mapping
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param alpha(in) : enhanced strain scalars
   * @return Core::LinAlg::Matrix<num_str, 1> : Enhanced Green-Lagrange strains E^{enh}
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
  evaluate_gl_strains_with_eas(
      const Discret::Elements::SpatialMaterialMapping<celltype>& displacement_based_mapping,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Discret::Elements::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const Core::LinAlg::Matrix<Discret::Elements::EasTypeToNumEas<eastype>::num_eas, 1>& alpha)
  {
    const Core::LinAlg::SymmetricTensor<double, Discret::Elements::num_dim<celltype>,
        Discret::Elements::num_dim<celltype>>
        displacement_based_cauchygreen =
            Discret::Elements::evaluate_cauchy_green<celltype>(displacement_based_mapping);

    Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        gl_strain =
            Discret::Elements::evaluate_green_lagrange_strain(displacement_based_cauchygreen) +
            evaluate_enhanced_assumed_gl_strains<celltype, eastype>(Mtilde, alpha);

    return gl_strain;
  }

  /*!
   * @brief Integrate the EAS stiffness matrices
   *
   * @tparam celltype, eastype
   * @param stress(in) : 2. Piola Kirchhoff stress tensor and material tangent
   * @param Mtilde(in) : matrix Mtilde in the material configuration
   * @param Bop(in) : B-operator
   * @param integration_factor(in) : integration factor (Gauss point weight times Jacobi
   * determinant)
   * @param eas_iteration_data(in/out) : EAS matrices and vectors
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  void integrate_eas(const Discret::Elements::Stress<celltype>& stress,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Discret::Elements::EasTypeToNumEas<eastype>::num_eas>& Mtilde,
      const Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
          Discret::Elements::num_dof_per_ele<celltype>>& Bop,
      const double integration_factor,
      Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data)
  {
    // integrate Kaa: Kaa += (Mtilde^T . cmat . Mtilde) * detJ * w(gp)
    // IMPORTANT: We save this in invKaa here since after the loop over all Gauss points, we
    // invert the matrix. At this point, this is still Kaa and NOT invKaa.
    Core::LinAlg::Matrix<Discret::Elements::num_str<celltype>,
        Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
        cmatM(Core::LinAlg::Initialization::zero);
    cmatM.multiply(Core::LinAlg::make_stress_like_voigt_view(stress.cmat_), Mtilde);
    eas_iteration_data.invKaa.multiply_tn(integration_factor, Mtilde, cmatM, 1.);

    // integrate Kda: Kda += (B^T . cmat . Mtilde) * detJ * w(gp)
    eas_iteration_data.Kda.multiply_tn(integration_factor, Bop, cmatM, 1.);

    // integrate s: s += (Mtilde^T . S) * detJ * w(gp)
    eas_iteration_data.s.multiply_tn(
        integration_factor, Mtilde, Core::LinAlg::make_stress_like_voigt_view(stress.pk2_), 1.);
  }

  /*!
   * @brief Add EAS internal force contribution of one Gauss point
   *
   * The EAS internal force contribution is $f_{eas} = - K_{da} K_{aa}^{-1} s$.
   *
   * @tparam celltype : Cell type
   * @param minusKdainvKaa(in) : matrix product $- K_{da} K_{aa}^{-1}$
   * @param s(in) : enhancement vector s
   * @param force(in/out) : internal force vector where the contribution is added to
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  void add_eas_internal_force(
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>,
          Discret::Elements::EasTypeToNumEas<eastype>::num_eas>& minusKdainvKaa,
      const Core::LinAlg::Matrix<Discret::Elements::EasTypeToNumEas<eastype>::num_eas, 1>& s,
      Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>, 1>& force_vector)
  {
    force_vector.multiply_nn(1.0, minusKdainvKaa, s, 1.0);
  }

  /*!
   * @brief Add EAS stiffness matrix contribution of one Gauss point
   *
   * The EAS stiffness matrix contribution is $- K_{da} K_{aa}^{-1} K_{da}^T$.
   *
   * @tparam celltype : Cell type
   * @param minusKdainvKaa(in) : matrix product $- K_{da} K_{aa}^{-1}$
   * @param Kda(in) : EAS stiffness matrix part K_{da}
   * @param stiffness_matrix(in/out) : stiffness matrix where the local contribution is added to
   */
  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  void add_eas_stiffness_matrix(
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>,
          Discret::Elements::EasTypeToNumEas<eastype>::num_eas>& minusKdainvKaa,
      const Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>,
          Discret::Elements::EasTypeToNumEas<eastype>::num_eas>& Kda,
      Core::LinAlg::Matrix<Discret::Elements::num_dof_per_ele<celltype>,
          Discret::Elements::num_dof_per_ele<celltype>>& stiffness_matrix)
  {
    stiffness_matrix.multiply_nt(1.0, minusKdainvKaa, Kda, 1.0);
  }

  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
  struct EASKinematics
  {
    static constexpr int num_str = Core::FE::dim<celltype> * (Core::FE::dim<celltype> + 1) / 2;
    static constexpr int num_dof_per_ele = Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>;
    Core::LinAlg::Matrix<num_str, num_dof_per_ele> b_op{};
    Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        enhanced_gl{};
    Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        enhanced_deformation_gradient{};
    Core::LinAlg::Matrix<num_str, Discret::Elements::EasTypeToNumEas<eastype>::num_eas> m_tilde{};
  };

  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
    requires(Core::FE::dim<celltype> == 2)
  struct AxisymmetricEASKinematics
  {
    AxisymmetricKinematics<celltype> displacement_based;
    Core::LinAlg::Tensor<double, 3, 3> enhanced_deformation_gradient;
    Core::LinAlg::SymmetricTensor<double, 3, 3> enhanced_gl;
    Core::LinAlg::Matrix<num_str<celltype>, EasTypeToNumEas<eastype>::num_eas> m_tilde;
  };

  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
    requires(Core::FE::dim<celltype> == 2)
  AxisymmetricEASKinematics<celltype, eastype> evaluate_axisymmetric_eas_kinematics(
      const ElementNodes<celltype>& nodal_coordinates,
      const CentroidTransformation<celltype>& centroid_transformation,
      const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const EasIterationData<celltype, eastype>& eas_iteration_data)
  {
    AxisymmetricEASKinematics<celltype, eastype> kinematics{
        .displacement_based = evaluate_axisymmetric_kinematics(
            nodal_coordinates, shape_functions, jacobian_mapping, false),
        .enhanced_deformation_gradient = {},
        .enhanced_gl = {},
        .m_tilde = evaluate_eas_shape_functions_material_config<celltype, eastype>(
            jacobian_mapping.determinant_, centroid_transformation, xi)};

    kinematics.enhanced_gl = kinematics.displacement_based.green_lagrange_strain;
    const auto enhanced_in_plane_strain = evaluate_enhanced_assumed_gl_strains<celltype, eastype>(
        kinematics.m_tilde, eas_iteration_data.alpha);
    for (int i = 0; i < 2; ++i)
      for (int j = i; j < 2; ++j) kinematics.enhanced_gl(i, j) += enhanced_in_plane_strain(i, j);

    kinematics.enhanced_deformation_gradient = compute_deformation_gradient_from_gl_strains(
        kinematics.displacement_based.deformation_gradient, kinematics.enhanced_gl);
    return kinematics;
  }

  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
    requires(Core::FE::dim<celltype> == 2)
  Core::LinAlg::SymmetricTensor<double, 3, 3> evaluate_axisymmetric_eas_strain_derivative(
      const Core::LinAlg::Matrix<num_str<celltype>, EasTypeToNumEas<eastype>::num_eas>& m_tilde,
      const int enhanced_strain_parameter)
  {
    Core::LinAlg::SymmetricTensor<double, 3, 3> derivative{};
    derivative(0, 0) = m_tilde(0, enhanced_strain_parameter);
    derivative(1, 1) = m_tilde(1, enhanced_strain_parameter);
    derivative(0, 1) = 0.5 * m_tilde(2, enhanced_strain_parameter);
    return derivative;
  }

  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype>
    requires(Core::FE::dim<celltype> == 2)
  void integrate_axisymmetric_eas(const AxisymmetricEASKinematics<celltype, eastype>& kinematics,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
      const double integration_factor, EasIterationData<celltype, eastype>& eas_iteration_data)
  {
    constexpr int num_dofs = num_dof_per_ele<celltype>;
    constexpr int num_eas = EasTypeToNumEas<eastype>::num_eas;
    std::array<Core::LinAlg::SymmetricTensor<double, 3, 3>, num_dofs>
        displacement_strain_derivatives;
    std::array<Core::LinAlg::SymmetricTensor<double, 3, 3>, num_eas> enhanced_strain_derivatives;

    for (int node = 0; node < num_nodes<celltype>; ++node)
      for (int component = 0; component < 2; ++component)
        displacement_strain_derivatives[2 * node + component] =
            evaluate_axisymmetric_strain_derivative(node, component, kinematics.displacement_based,
                shape_functions, jacobian_mapping, false);
    for (int parameter = 0; parameter < num_eas; ++parameter)
      enhanced_strain_derivatives[parameter] =
          evaluate_axisymmetric_eas_strain_derivative<celltype, eastype>(
              kinematics.m_tilde, parameter);

    auto contract_with_material_tangent =
        [&](const Core::LinAlg::SymmetricTensor<double, 3, 3>& left,
            const Core::LinAlg::SymmetricTensor<double, 3, 3>& right)
    {
      double value = 0.0;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < 3; ++k)
            for (int l = 0; l < 3; ++l)
              value += left(i, j) * stress.cmat_3d_(i, j, k, l) * right(k, l);
      return value;
    };

    for (int row = 0; row < num_eas; ++row)
    {
      double enhanced_residual = 0.0;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          enhanced_residual += enhanced_strain_derivatives[row](i, j) * stress.pk2_3d_(i, j);
      eas_iteration_data.s(row) += integration_factor * enhanced_residual;

      for (int column = 0; column < num_eas; ++column)
        eas_iteration_data.invKaa(row, column) +=
            integration_factor * contract_with_material_tangent(enhanced_strain_derivatives[row],
                                     enhanced_strain_derivatives[column]);

      for (int displacement_dof = 0; displacement_dof < num_dofs; ++displacement_dof)
        eas_iteration_data.Kda(displacement_dof, row) +=
            integration_factor *
            contract_with_material_tangent(displacement_strain_derivatives[displacement_dof],
                enhanced_strain_derivatives[row]);
    }
  }

  template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
      FourC::Solid::KinemType kinematic_type>
  EASKinematics<celltype, eastype> evaluate_eas_kinematics(
      const Discret::Elements::ElementNodes<celltype> nodal_coordinates,
      const Discret::Elements::CentroidTransformation<celltype>& centeroid_transformation,
      const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
      const Discret::Elements::JacobianMapping<celltype>& jacobian_mapping,
      const Discret::Elements::EasIterationData<celltype, eastype>& eas_iteration_data)
  {
    EASKinematics<celltype, eastype> eas_kinematics{};

    if constexpr (kinematic_type == FourC::Solid::KinemType::nonlinearTotLag)
    {
      const Discret::Elements::SpatialMaterialMapping<celltype>
          displacement_based_spatial_material_mapping =
              evaluate_spatial_material_mapping(jacobian_mapping, nodal_coordinates);

      eas_kinematics.b_op = Discret::Elements::evaluate_strain_gradient(
          jacobian_mapping, displacement_based_spatial_material_mapping);

      eas_kinematics.m_tilde = evaluate_eas_shape_functions_material_config<celltype, eastype>(
          jacobian_mapping.determinant_, centeroid_transformation, xi);

      eas_kinematics.enhanced_gl = evaluate_gl_strains_with_eas<celltype, eastype>(
          displacement_based_spatial_material_mapping, eas_kinematics.m_tilde,
          eas_iteration_data.alpha);

      eas_kinematics.enhanced_deformation_gradient =
          Discret::Elements::compute_deformation_gradient_from_gl_strains(
              displacement_based_spatial_material_mapping.deformation_gradient_,
              eas_kinematics.enhanced_gl);
    }
    else if constexpr (kinematic_type == FourC::Solid::KinemType::linear)
    {
      eas_kinematics.b_op = Discret::Elements::evaluate_linear_strain_gradient(jacobian_mapping);

      eas_kinematics.m_tilde = evaluate_eas_shape_functions_material_config<celltype, eastype>(
          jacobian_mapping.determinant_, centeroid_transformation, xi);

      Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
          gl_strain_displacement_based =
              evaluate_linear_gl_strain(nodal_coordinates, eas_kinematics.b_op);

      eas_kinematics.enhanced_gl = gl_strain_displacement_based;
      eas_kinematics.enhanced_gl += evaluate_enhanced_assumed_gl_strains<celltype, eastype>(
          eas_kinematics.m_tilde, eas_iteration_data.alpha);

      eas_kinematics.enhanced_deformation_gradient =
          Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double,
              Core::FE::dim<celltype>, Core::FE::dim<celltype>>);
    }
    else
    {
      FOUR_C_THROW("Unknown kinematic type!");
    }

    return eas_kinematics;
  }
}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif
