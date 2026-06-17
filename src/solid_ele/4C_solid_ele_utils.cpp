// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solid_ele_utils.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io_input_parameter_container.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_symmetric_tensor_eigen.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_solid_ele_calc_lib.hpp"
#include "4C_solid_ele_properties.hpp"

#include <algorithm>

FOUR_C_NAMESPACE_OPEN



int Solid::Utils::ReadElement::read_element_material(
    const Core::IO::InputParameterContainer& container)
{
  int material = container.get<int>("MAT");
  return material;
}


template <unsigned dim>
Discret::Elements::SolidElementProperties<dim>
Solid::Utils::ReadElement::read_solid_element_properties(
    const Core::IO::InputParameterContainer& container)
{
  Discret::Elements::SolidElementProperties<dim> solid_properties{};

  // element technology
  solid_properties.element_technology = container.get_or<Discret::Elements::ElementTechnology>(
      "TECH", Discret::Elements::ElementTechnology::none);

  // prestress technology
  solid_properties.prestress_technology = container.get_or<Discret::Elements::PrestressTechnology>(
      "PRESTRESS_TECH", Discret::Elements::PrestressTechnology::none);

  // kinematic type
  solid_properties.kintype = container.get<Inpar::Solid::KinemType>("KINEM");

  if constexpr (dim == 2)
  {
    solid_properties.reference_thickness = container.get_or<double>("THICKNESS", 1.0);
    solid_properties.plane_assumption =
        container.get<Discret::Elements::PlaneAssumption>("PLANE_ASSUMPTION");
  }

  return solid_properties;
}

void Solid::Utils::nodal_block_information_solid(
    Core::Elements::Element* dwele, int& numdf, int& dimns)
{
  numdf = Core::FE::get_dimension(dwele->shape());
  dimns = numdf * (numdf + 1) / 2;
}


template Discret::Elements::SolidElementProperties<2>
Solid::Utils::ReadElement::read_solid_element_properties(
    const Core::IO::InputParameterContainer& container);
template Discret::Elements::SolidElementProperties<3>
Solid::Utils::ReadElement::read_solid_element_properties(
    const Core::IO::InputParameterContainer& container);

FOUR_C_NAMESPACE_CLOSE
