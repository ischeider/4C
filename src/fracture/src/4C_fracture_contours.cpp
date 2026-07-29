// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fracture_contours.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

Fracture::ContourLayers Fracture::build_contour_layers(
    const std::vector<VolumeElementTopology>& elements, const std::set<int>& crack_front_node_ids,
    const int number_of_contours)
{
  TEUCHOS_FUNC_TIME_MONITOR("Fracture: find domains");

  FOUR_C_ASSERT_ALWAYS(number_of_contours > 0, "The number of fracture contours must be positive.");
  FOUR_C_ASSERT_ALWAYS(
      !crack_front_node_ids.empty(), "Cannot build contours for an empty crack front.");

  std::map<int, const VolumeElementTopology*> elements_by_id;
  std::map<std::vector<int>, std::vector<int>> elements_by_face;
  std::map<int, std::vector<int>> elements_by_node;
  ContourLayers result;
  for (const auto& element : elements)
  {
    FOUR_C_ASSERT_ALWAYS(elements_by_id.emplace(element.element_id, &element).second,
        "Duplicate element ID {} in fracture contour topology.", element.element_id);
    result.element_layer[element.element_id] = 0;
    FOUR_C_ASSERT_ALWAYS(
        !element.node_ids.empty(), "Element {} contains no nodes.", element.element_id);
    auto sorted_node_ids = element.node_ids;
    std::ranges::sort(sorted_node_ids);
    sorted_node_ids.erase(std::ranges::unique(sorted_node_ids).begin(), sorted_node_ids.end());
    for (const int node_id : sorted_node_ids)
      elements_by_node[node_id].push_back(element.element_id);
    for (auto face : element.faces)
    {
      FOUR_C_ASSERT_ALWAYS(!face.empty(), "Element {} contains an empty face.", element.element_id);
      std::ranges::sort(face);
      face.erase(std::ranges::unique(face).begin(), face.end());
      if (face.size() < 3) continue;
      elements_by_face[std::move(face)].push_back(element.element_id);
    }
  }

  std::map<int, std::set<int>> neighbors;
  for (const auto& [node, incident_elements] : elements_by_node)
    for (const int element_id : incident_elements)
      neighbors[element_id].insert(incident_elements.begin(), incident_elements.end());
  for (auto& [element_id, element_neighbors] : neighbors) element_neighbors.erase(element_id);

  for (const auto& [face, incident_elements] : elements_by_face)
  {
    FOUR_C_ASSERT_ALWAYS(incident_elements.size() <= 2,
        "Non-manifold face with {} incident elements in fracture contour topology.",
        incident_elements.size());
  }

  std::set<int> current_layer;
  for (const auto& element : elements)
  {
    if (std::ranges::any_of(element.node_ids,
            [&](const int node_id) { return crack_front_node_ids.contains(node_id); }))
    {
      current_layer.insert(element.element_id);
      result.element_layer[element.element_id] = 1;
    }
  }
  FOUR_C_ASSERT_ALWAYS(
      !current_layer.empty(), "No volume element is attached to the configured crack front.");

  std::set<int> cumulative = current_layer;
  result.contours.reserve(number_of_contours);
  result.contours.push_back(cumulative);

  for (int layer = 2; layer <= number_of_contours; ++layer)
  {
    std::set<int> next_layer;
    for (const int element_id : current_layer)
    {
      for (const int neighbor_id : neighbors[element_id])
      {
        if (result.element_layer.at(neighbor_id) != 0) continue;
        result.element_layer[neighbor_id] = layer;
        next_layer.insert(neighbor_id);
      }
    }
    cumulative.insert(next_layer.begin(), next_layer.end());
    result.contours.push_back(cumulative);
    current_layer = std::move(next_layer);
  }

  return result;
}

void Fracture::normalize_extension_fields(std::vector<std::map<int, double>>& extension_fields)
{
  FOUR_C_ASSERT_ALWAYS(
      !extension_fields.empty(), "Cannot normalize an empty set of fracture extension fields.");

  for (const auto& [node, value] : extension_fields.front())
  {
    double partition_sum = 0.0;
    for (const auto& field : extension_fields)
    {
      FOUR_C_ASSERT_ALWAYS(field.size() == extension_fields.front().size(),
          "Fracture extension fields are defined on different node sets.");
      const auto entry = field.find(node);
      FOUR_C_ASSERT_ALWAYS(
          entry != field.end(), "Fracture extension fields are defined on different node sets.");
      partition_sum += entry->second;
    }
    if (partition_sum <= 1.0e-14) continue;
    for (auto& field : extension_fields) field.at(node) /= partition_sum;
  }
}

FOUR_C_NAMESPACE_CLOSE
