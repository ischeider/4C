// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_CALC_LIB_PLANE_HPP
#define FOUR_C_SOLID_ELE_CALC_LIB_PLANE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_solid_ele_calc_lib.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  void transform_to_3d(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, int gp, int eleGID,
      const std::function<void(const Core::LinAlg::Tensor<double, 3, 3>&,
          const Core::LinAlg::SymmetricTensor<double, 3, 3>&, const Mat::EvaluationContext<3>&)>&
          funct,
      const OutOfPlaneKinematics& out_of_plane = {});

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  Stress<celltype> evaluate_material_stress(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID,
      const OutOfPlaneKinematics& out_of_plane = {});

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  void update_material(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID,
      const OutOfPlaneKinematics& out_of_plane = {});

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  [[nodiscard]] double evaluate_material_strain_energy(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID,
      const OutOfPlaneKinematics& out_of_plane = {});

  /*!
   * @brief Computes the out-of-plane (hoop) kinematics for the axisymmetric formulation.
   *
   * The hoop stretch is F_theta = r_current / r_reference = 1 + u_r / r_reference. For linear
   * kinematics the hoop Green-Lagrange strain reduces to u_r / r_reference, while for nonlinear
   * kinematics it is the consistent finite strain 0.5 * (F_theta^2 - 1).
   *
   * @param kintype (in) : Kinematic type (linear or nonlinear)
   * @param reference_radius (in) : Radial reference coordinate r of the Gauss point (> 0)
   * @param radial_displacement (in) : Radial displacement u_r at the Gauss point
   */
  inline OutOfPlaneKinematics compute_axisymmetric_out_of_plane(
      Inpar::Solid::KinemType kintype, double reference_radius, double radial_displacement)
  {
    const double hoop_stretch = 1.0 + radial_displacement / reference_radius;
    if (kintype == Inpar::Solid::KinemType::linear)
    {
      return {.gl_strain_33 = radial_displacement / reference_radius, .defgrd_33 = 1.0};
    }
    return {.gl_strain_33 = 0.5 * (hoop_stretch * hoop_stretch - 1.0), .defgrd_33 = hoop_stretch};
  }

  /*!
   * @brief Adds the internal force and stiffness contributions of the out-of-plane (hoop)
   * direction for the axisymmetric formulation at one Gauss point.
   *
   * The in-plane contributions are handled by the regular element formulation. This function adds
   * the additional terms that arise because the hoop strain E_theta couples to the radial
   * displacement and, through the material, to the in-plane stresses.
   *
   * @tparam celltype : Cell type (2D)
   * @param shape_functions (in) : Shape functions evaluated at the Gauss point
   * @param jacobian_mapping (in) : Jacobian mapping (provides spatial shape function derivatives)
   * @param defgrd (in) : In-plane deformation gradient
   * @param stress (in) : Stress measures including the full 3D PK2 stress and material tangent
   * @param reference_radius (in) : Radial reference coordinate r of the Gauss point (> 0)
   * @param hoop_stretch (in) : Out-of-plane stretch F_theta
   * @param kintype (in) : Kinematic type (linear or nonlinear)
   * @param integration_fac (in) : Integration factor (already includes the r weight)
   * @param force_vector (in/out) : Optional element force vector
   * @param stiffness_matrix (in/out) : Optional element stiffness matrix
   */
  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  void add_axisymmetric_hoop_force_stiffness(
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const JacobianMapping<celltype>& jacobian_mapping,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd, const Stress<celltype>& stress,
      const double reference_radius, const double hoop_stretch, Inpar::Solid::KinemType kintype,
      const double integration_fac,
      Core::LinAlg::Matrix<2 * Core::FE::num_nodes(celltype), 1>* force_vector,
      Core::LinAlg::Matrix<2 * Core::FE::num_nodes(celltype), 2 * Core::FE::num_nodes(celltype)>*
          stiffness_matrix)
  {
    constexpr std::size_t num_nodes = Core::FE::num_nodes(celltype);
    const bool is_linear = (kintype == Inpar::Solid::KinemType::linear);

    // Out-of-plane PK2 stress S_theta and material tangents (full 3D, index 2 = hoop).
    const double s_theta = stress.pk2_3d_(2, 2);
    const double c_theta_theta = stress.cmat_3d_(2, 2, 2, 2);

    // In-plane deformation gradient (identity for the linear case).
    const Core::LinAlg::Tensor<double, 2, 2> F =
        is_linear ? Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 2, 2>)
                  : defgrd;

    // Coupling of the hoop direction to the in-plane material tangent: C_theta_inplane(a,b) =
    // cmat_3d_(a,b,2,2).
    auto compute_w = [&](std::size_t node) -> std::array<double, 2>
    {
      // c_grad(a) = sum_b cmat_3d_(a,b,theta,theta) * dN_i/dX_b
      std::array<double, 2> c_grad{};
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          c_grad[a] += stress.cmat_3d_(a, b, 2, 2) * jacobian_mapping.N_XYZ[node](b);

      // w(k) = sum_a F(k,a) * c_grad(a)
      std::array<double, 2> w{};
      for (int k = 0; k < 2; ++k)
        for (int a = 0; a < 2; ++a) w[k] += F(k, a) * c_grad[a];
      return w;
    };

    for (std::size_t i = 0; i < num_nodes; ++i)
    {
      const double N_i = shape_functions.shapefunctions_(i);
      // Hoop strain-displacement operator (radial dof only): B_theta_i = F_theta * N_i / r.
      const double B_theta_i = hoop_stretch * N_i / reference_radius;

      if (force_vector != nullptr)
      {
        (*force_vector)(2 * i + 0) += s_theta * B_theta_i * integration_fac;
      }

      if (stiffness_matrix != nullptr)
      {
        const std::array<double, 2> w_i = compute_w(i);
        for (std::size_t j = 0; j < num_nodes; ++j)
        {
          const double N_j = shape_functions.shapefunctions_(j);
          const double B_theta_j = hoop_stretch * N_j / reference_radius;

          // Hoop-hoop material stiffness.
          (*stiffness_matrix)(2 * i + 0, 2 * j + 0) +=
              B_theta_i * c_theta_theta * B_theta_j * integration_fac;

          // Geometric hoop stiffness (nonlinear only): d(delta E_theta) term.
          if (!is_linear)
          {
            (*stiffness_matrix)(2 * i + 0, 2 * j + 0) +=
                s_theta * (N_i / reference_radius) * (N_j / reference_radius) * integration_fac;
          }

          // Coupling between in-plane displacements and the hoop stress, and its transpose.
          for (int k = 0; k < 2; ++k)
          {
            (*stiffness_matrix)(2 * i + k, 2 * j + 0) += w_i[k] * B_theta_j * integration_fac;
            (*stiffness_matrix)(2 * j + 0, 2 * i + k) += B_theta_j * w_i[k] * integration_fac;
          }
        }
      }
    }
  }

}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif