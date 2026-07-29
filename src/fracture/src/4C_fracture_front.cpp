// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fracture_front.hpp"

#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

FOUR_C_NAMESPACE_OPEN

namespace
{
  using Fracture::FrontSegment;
  using Fracture::Vector3;

  Vector3 add(const Vector3& left, const Vector3& right)
  {
    return {left[0] + right[0], left[1] + right[1], left[2] + right[2]};
  }

  Vector3 scale(const Vector3& vector, const double factor)
  {
    return {factor * vector[0], factor * vector[1], factor * vector[2]};
  }

  double norm(const Vector3& vector)
  {
    return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
  }

  Vector3 normalize(const Vector3& vector)
  {
    const double length = norm(vector);
    FOUR_C_ASSERT_ALWAYS(length > 100.0 * std::numeric_limits<double>::epsilon(),
        "Degenerate crack-front segment or tangent.");
    return scale(vector, 1.0 / length);
  }

  std::vector<double> shape_functions(const std::size_t number_of_nodes, const double xi)
  {
    if (number_of_nodes == 2) return {0.5 * (1.0 - xi), 0.5 * (1.0 + xi)};
    if (number_of_nodes == 3) return {0.5 * xi * (xi - 1.0), 1.0 - xi * xi, 0.5 * xi * (xi + 1.0)};
    FOUR_C_THROW("Crack-front segments support only line2 and line3 interpolation.");
  }

  std::vector<double> shape_derivatives(const std::size_t number_of_nodes, const double xi)
  {
    if (number_of_nodes == 2) return {-0.5, 0.5};
    if (number_of_nodes == 3) return {xi - 0.5, -2.0 * xi, xi + 0.5};
    FOUR_C_THROW("Crack-front segments support only line2 and line3 interpolation.");
  }

  Vector3 derivative(
      const FrontSegment& segment, const std::map<int, Vector3>& coordinates, const double xi)
  {
    const auto derivatives = shape_derivatives(segment.node_ids.size(), xi);
    Vector3 result{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < segment.node_ids.size(); ++i)
      result = add(result, scale(coordinates.at(segment.node_ids[i]), derivatives[i]));
    return result;
  }

  void integrate_support(const FrontSegment& segment, const std::map<int, Vector3>& coordinates,
      std::map<int, double>& support)
  {
    constexpr std::array<double, 5> points = {
        -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
    constexpr std::array<double, 5> weights = {0.2369268850561891, 0.4786286704993665,
        0.5688888888888889, 0.4786286704993665, 0.2369268850561891};

    for (std::size_t gp = 0; gp < points.size(); ++gp)
    {
      const auto shapes = shape_functions(segment.node_ids.size(), points[gp]);
      const double jacobian = norm(derivative(segment, coordinates, points[gp]));
      FOUR_C_ASSERT_ALWAYS(jacobian > 100.0 * std::numeric_limits<double>::epsilon(),
          "Degenerate crack-front segment.");
      for (std::size_t i = 0; i < segment.node_ids.size(); ++i)
        support[segment.node_ids[i]] += weights[gp] * shapes[i] * jacobian;
    }
  }

  std::vector<double> nodal_coordinates(const std::size_t number_of_nodes)
  {
    if (number_of_nodes == 2) return {-1.0, 1.0};
    if (number_of_nodes == 3) return {-1.0, 0.0, 1.0};
    FOUR_C_THROW("Crack-front segments support only line2 and line3 interpolation.");
  }
}  // namespace

Fracture::Front Fracture::build_front(
    const std::vector<FrontSegment>& segments, const std::map<int, Vector3>& coordinates)
{
  FOUR_C_ASSERT_ALWAYS(!segments.empty(), "A crack front requires at least one line segment.");

  std::map<int, std::vector<std::size_t>> endpoint_segments;
  for (std::size_t segment_index = 0; segment_index < segments.size(); ++segment_index)
  {
    const auto& node_ids = segments[segment_index].node_ids;
    FOUR_C_ASSERT_ALWAYS(node_ids.size() == 2 || node_ids.size() == 3,
        "Crack-front segment {} must be line2 or line3.", segment_index);
    FOUR_C_ASSERT_ALWAYS(node_ids.front() != node_ids.back(),
        "Crack-front segment {} has identical endpoints.", segment_index);
    for (const int node_id : node_ids)
      FOUR_C_ASSERT_ALWAYS(
          coordinates.contains(node_id), "Missing coordinates for crack-front node {}.", node_id);
    endpoint_segments[node_ids.front()].push_back(segment_index);
    endpoint_segments[node_ids.back()].push_back(segment_index);
  }

  std::vector<int> endpoints;
  for (const auto& [node_id, incident_segments] : endpoint_segments)
  {
    FOUR_C_ASSERT_ALWAYS(incident_segments.size() <= 2,
        "Crack-front node {} branches into {} segments.", node_id, incident_segments.size());
    if (incident_segments.size() == 1) endpoints.push_back(node_id);
  }
  FOUR_C_ASSERT_ALWAYS(endpoints.empty() || endpoints.size() == 2,
      "A crack front must be one open path or one closed cycle; found {} endpoints.",
      endpoints.size());

  const bool closed = endpoints.empty();
  int current_node = closed ? endpoint_segments.begin()->first
                            : *std::min_element(endpoints.begin(), endpoints.end());
  std::vector<bool> used(segments.size(), false);
  std::vector<FrontSegment> ordered_segments;
  ordered_segments.reserve(segments.size());
  std::vector<int> ordered_node_ids{current_node};
  std::set<int> visited_nodes{current_node};

  for (std::size_t step = 0; step < segments.size(); ++step)
  {
    const auto& candidates = endpoint_segments.at(current_node);
    std::optional<std::size_t> selected;
    int selected_other_node = std::numeric_limits<int>::max();
    for (const std::size_t candidate : candidates)
    {
      if (used[candidate]) continue;
      const auto& ids = segments[candidate].node_ids;
      const int other = ids.front() == current_node ? ids.back() : ids.front();
      if (!selected || other < selected_other_node)
      {
        selected = candidate;
        selected_other_node = other;
      }
    }
    FOUR_C_ASSERT_ALWAYS(selected.has_value(),
        "Crack-front line segments are disconnected or contain a premature cycle.");

    used[*selected] = true;
    FrontSegment segment = segments[*selected];
    if (segment.node_ids.front() != current_node) std::ranges::reverse(segment.node_ids);

    for (std::size_t i = 1; i < segment.node_ids.size(); ++i)
    {
      const int node_id = segment.node_ids[i];
      const bool closes_cycle = closed && step + 1 == segments.size() &&
                                i + 1 == segment.node_ids.size() &&
                                node_id == ordered_node_ids.front();
      FOUR_C_ASSERT_ALWAYS(closes_cycle || visited_nodes.insert(node_id).second,
          "Crack-front node {} occurs more than once in the ordered front.", node_id);
      if (!closes_cycle) ordered_node_ids.push_back(node_id);
    }
    current_node = segment.node_ids.back();
    ordered_segments.push_back(std::move(segment));
  }

  FOUR_C_ASSERT_ALWAYS(std::ranges::all_of(used, [](const bool value) { return value; }),
      "Crack-front line segments are disconnected.");
  FOUR_C_ASSERT_ALWAYS((closed && current_node == ordered_node_ids.front()) ||
                           (!closed && current_node != ordered_node_ids.front()),
      "Invalid crack-front closure.");

  std::map<int, double> support;
  std::map<int, Vector3> tangent_sums;
  for (const FrontSegment& segment : ordered_segments)
  {
    integrate_support(segment, coordinates, support);
    const auto xis = nodal_coordinates(segment.node_ids.size());
    for (std::size_t i = 0; i < segment.node_ids.size(); ++i)
      tangent_sums[segment.node_ids[i]] = add(
          tangent_sums[segment.node_ids[i]], normalize(derivative(segment, coordinates, xis[i])));
  }

  std::vector<FrontNode> nodes;
  nodes.reserve(ordered_node_ids.size());
  for (const int node_id : ordered_node_ids)
  {
    FOUR_C_ASSERT_ALWAYS(support.at(node_id) > 0.0,
        "Crack-front node {} has a nonpositive support measure.", node_id);
    nodes.push_back({.node_id = node_id,
        .tangent = normalize(tangent_sums.at(node_id)),
        .support_measure = support.at(node_id)});
  }

  return {.nodes = std::move(nodes), .segments = std::move(ordered_segments), .closed = closed};
}

FOUR_C_NAMESPACE_CLOSE
