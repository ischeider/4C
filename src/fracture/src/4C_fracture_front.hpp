// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FRACTURE_FRONT_HPP
#define FOUR_C_FRACTURE_FRONT_HPP

#include "4C_config.hpp"

#include <array>
#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Fracture
{
  using Vector3 = std::array<double, 3>;

  struct FrontSegment
  {
    /// Ordered line2 {end, end} or line3 {end, midpoint, end} node IDs.
    std::vector<int> node_ids;
  };

  struct FrontNode
  {
    int node_id;
    Vector3 tangent;
    double support_measure;
  };

  struct Front
  {
    std::vector<FrontNode> nodes;
    std::vector<FrontSegment> segments;
    bool closed;
  };

  /**
   * Order and validate one connected crack-front path or cycle.
   *
   * Support measures are evaluated as integral N_a ds on the supplied line interpolation using
   * Gaussian quadrature. Segment orientation is made consistent with the returned node ordering.
   */
  Front build_front(
      const std::vector<FrontSegment>& segments, const std::map<int, Vector3>& coordinates);
}  // namespace Fracture

FOUR_C_NAMESPACE_CLOSE

#endif
