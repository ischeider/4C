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
          funct);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  Stress<celltype> evaluate_material_stress(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  Stress<celltype> evaluate_axisymmetric_material_stress(Mat::So3Material& material,
      const Core::LinAlg::Tensor<double, 3, 3>& defgrd,
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, int gp, int eleGID);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  void update_material(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  void update_axisymmetric_material(Mat::So3Material& material,
      const Core::LinAlg::Tensor<double, 3, 3>& defgrd, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, int gp, int eleGID);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  [[nodiscard]] double evaluate_material_strain_energy(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  [[nodiscard]] double evaluate_axisymmetric_material_strain_energy(Mat::So3Material& material,
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& gl_strain,
      const Mat::EvaluationContext<2>& context, int gp, int eleGID);

}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif