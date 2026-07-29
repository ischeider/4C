// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fracture_crack_front.hpp"

#include "4C_fem_condition.hpp"
#include "4C_fracture_input.hpp"
#include "4C_io_input_parameter_container.hpp"
#include "4C_utils_exceptions.hpp"

#include <charconv>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace
{
  std::optional<int> parse_legacy_line_id(const std::string& crack_id)
  {
    int id = 0;
    const auto [end, error] =
        std::from_chars(crack_id.data(), crack_id.data() + crack_id.size(), id);
    if (error != std::errc{} || end != crack_id.data() + crack_id.size()) return std::nullopt;
    if (id <= 0) FOUR_C_THROW("Legacy fracture DLINE ID '{}' must be positive.", crack_id);
    return id - 1;
  }

  Core::IO::InputParameterContainer condition_data(
      const Fracture::Crack& crack, const std::size_t configuration_index)
  {
    Core::IO::InputParameterContainer data;
    data.add("configuration_index", static_cast<int>(configuration_index));
    data.add("crack_id", crack.crack_id);
    return data;
  }
}  // namespace

std::vector<std::unique_ptr<Core::Conditions::Condition>> Fracture::make_crack_front_conditions(
    const std::vector<Crack>& cracks, const Core::Conditions::InputNodeSets& input_node_sets)
{
  std::vector<std::unique_ptr<Core::Conditions::Condition>> conditions;
  conditions.reserve(cracks.size());

  for (std::size_t i = 0; i < cracks.size(); ++i)
  {
    const Crack& crack = cracks[i];
    if (const auto named_set = input_node_sets.node_sets_names.find(crack.crack_id);
        named_set != input_node_sets.node_sets_names.end())
    {
      FOUR_C_ASSERT_ALWAYS(named_set->second.size() == 1,
          "Fracture crack_id '{}' is not unique in the mesh file ({} occurrences).", crack.crack_id,
          named_set->second.size());
      const int node_set_id = named_set->second.front();
      FOUR_C_ASSERT_ALWAYS(input_node_sets.node_sets.contains(node_set_id),
          "Fracture crack_id '{}' refers to missing node set ID {}.", crack.crack_id, node_set_id);
      const auto& nodes = input_node_sets.node_sets.at(node_set_id);
      FOUR_C_ASSERT_ALWAYS(
          !nodes.empty(), "Fracture crack_id '{}' resolves to an empty node set.", crack.crack_id);

      conditions.push_back(
          std::make_unique<Core::Conditions::Condition>(node_set_id, Core::Conditions::none, false,
              Core::Conditions::geometry_type_line, Core::Conditions::EntityType::node_set_name,
              nodes, condition_data(crack, i), crack.crack_id));
      continue;
    }

    const std::optional<int> legacy_id = parse_legacy_line_id(crack.crack_id);
    FOUR_C_ASSERT_ALWAYS(legacy_id.has_value(),
        "Fracture crack_id '{}' is neither an external node-set name nor a legacy DLINE ID.",
        crack.crack_id);
    FOUR_C_ASSERT_ALWAYS(static_cast<std::size_t>(*legacy_id) < input_node_sets.dline_fenode.size(),
        "Fracture crack_id '{}' refers to missing legacy DLINE {}.", crack.crack_id,
        *legacy_id + 1);
    const auto& nodes = input_node_sets.dline_fenode[*legacy_id];
    FOUR_C_ASSERT_ALWAYS(
        !nodes.empty(), "Fracture crack_id '{}' resolves to an empty DLINE.", crack.crack_id);

    conditions.push_back(std::make_unique<Core::Conditions::Condition>(*legacy_id,
        Core::Conditions::none, false, Core::Conditions::geometry_type_line,
        Core::Conditions::EntityType::legacy_id, nodes, condition_data(crack, i)));
  }

  return conditions;
}

FOUR_C_NAMESPACE_CLOSE
