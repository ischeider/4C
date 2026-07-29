// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FRACTURE_CONTOURS_HPP
#define FOUR_C_FRACTURE_CONTOURS_HPP

#include "4C_config.hpp"

#include <map>
#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Fracture
{
  struct VolumeElementTopology
  {
    int element_id;
    std::vector<int> node_ids;
    std::vector<std::vector<int>> faces;
  };

  struct ContourLayers
  {
    /// First topological layer for every supplied element; zero means outside requested contours.
    std::map<int, int> element_layer;

    /// Cumulative element IDs in contours 1..number_of_contours.
    std::vector<std::set<int>> contours;
  };

  /**
   * Build topological EDI contours using shared-node element adjacency.
   *
   * Layer one contains every element incident to a crack-front node. Each later layer contains
   * every unassigned element sharing at least one node with the preceding layer. No assumption is
   * made about closed rings, so symmetry and free-surface truncations require no special handling.
   */
  ContourLayers build_contour_layers(const std::vector<VolumeElementTopology>& elements,
      const std::set<int>& crack_front_node_ids, int number_of_contours);

  /**
   * Normalize crack-front extension weights to a partition of unity at every interior node.
   *
   * Nodes where all fields are zero remain zero and form the artificial contour boundary.
   */
  void normalize_extension_fields(std::vector<std::map<int, double>>& extension_fields);
}  // namespace Fracture

FOUR_C_NAMESPACE_CLOSE

#endif
