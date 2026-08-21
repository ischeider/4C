// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_NEUMANN_EVALUATOR_HPP
#define FOUR_C_SOLID_ELE_NEUMANN_EVALUATOR_HPP
#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret::Elements
{
  /*!
   * @brief Evaluates a Neumann condition @p condition for the element @p element.
   *
   * The element force vector is
   *
   * \f[
   * \boldsymbol{f}^{(e)} = \left[
   *    f_x^{1(e)}~f_y^{1(e)}~f_z^{1(e)}~\cdots~f_x^{n(e)}~f_y^{n(e)}~f_z^{n(e)}
   * \right]
   * @f]
   * with
   * @f[
   *   f_{x/y/z}^{i(e)} = \int_{\Omega^{(e)}} N^i \cdot \mathrm{value}_{x/y/z} \cdot
   *   \mathrm{funct}_{x/y/z} (t)  \mathrm{d} \Omega
   * \f],
   * where \f$n\f$ is the number of nodes of the element and \f$N^i\f$ is the \f$i\f$-th shape
   * function of the element.
   *
   * @note This function determines the shape of the element at runtime and calls the respective
   * templated version of @p evaluate_neumann. If you already know the Core::FE::CellType
   * of the element at compile-time, you could directly call @evaluate_neumann.
   *
   * @param element (in) : The element where we integrate
   * @param discretization (in) : discretization
   * @param condition (in) : The Neumann condition to be evaluated within the element.
   * @param dof_index_array (in) : The index array of the dofs of the element
   * @param element_force_vector (out) : The element force vector for the evaluated Neumann
   * condition
   * @param total_time (in) : The total time for time dependent Neumann conditions
   */
  template <int dim>
  void evaluate_neumann_by_element(Core::Elements::Element& element,
      const Core::FE::Discretization& discretization, const Core::Conditions::Condition& condition,
      Core::LinAlg::SerialDenseVector& element_force_vector, double total_time,
      bool axisymmetric = false);

  /*!
   * @brief Evaluates a Neumann condition @p condition for the element @p element with the
   * cell type known at compile time.
   *
   * The element force vector is
   *
   * @f[
   * \boldsymbol{f}_{(e)} = \left[
   *    f_x^{1(e)}~f_y^{1(e)}~f_z^{1(e)}~\cdots~f_x^{n(e)}~f_y^{n(e)}~f_z^{n(e)}
   * \right]
   * @f]
   * with
   * @f[
   * f_{x/y/z}^{i(e)} = \int_{\Omega^{(e)}} N^i \cdot \mathrm{value}_{x/y/z} \cdot
   * \mathrm{funct}_{x/y/z} (t)  \mathrm{d} \Omega
   * @f],
   * where \f$n\f$ is the number of nodes of the element and \f$N^i\f$ is the \f$i\f$-th shape
   * function of the element.
   *
   * @tparam celltype Cell type known at compile time
   *
   * @param element (in) : The element where we integrate
   * @param discretization (in) : discretization
   * @param condition (in) : The Neumann condition to be evaluated within the element.
   * @param dof_index_array (in) : The index array of the dofs of the element
   * @param element_force_vector (out) : The element force vector for the evaluated Neumann
   * condition
   * @param total_time (in) : The total time for time dependent Neumann conditions
   */
  template <Core::FE::CellType celltype, int dim>
  void evaluate_neumann(Core::Elements::Element& element,
      const Core::FE::Discretization& discretization, const Core::Conditions::Condition& condition,
      Core::LinAlg::SerialDenseVector& element_force_vector, double total_time,
      bool axisymmetric = false);

}  // namespace Discret::Elements


FOUR_C_NAMESPACE_CLOSE

#endif
