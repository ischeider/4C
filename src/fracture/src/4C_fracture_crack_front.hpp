// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FRACTURE_CRACK_FRONT_HPP
#define FOUR_C_FRACTURE_CRACK_FRONT_HPP

#include "4C_config.hpp"

#include <memory>
#include <string_view>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Conditions
{
  class Condition;
  struct InputNodeSets;
}  // namespace Core::Conditions

namespace Fracture
{
  struct Crack;

  inline constexpr std::string_view crack_front_condition_name = "CrackFront";

  /**
   * Resolve crack identifiers against external named node sets or legacy DLINE entities.
   *
   * External names take precedence. A decimal crack identifier is interpreted as a one-based
   * legacy DLINE ID only when no external node set with that exact name exists.
   */
  std::vector<std::unique_ptr<Core::Conditions::Condition>> make_crack_front_conditions(
      const std::vector<Crack>& cracks, const Core::Conditions::InputNodeSets& input_node_sets);
}  // namespace Fracture

FOUR_C_NAMESPACE_CLOSE

#endif
