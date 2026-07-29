// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FRACTURE_PROPAGATION_HPP
#define FOUR_C_FRACTURE_PROPAGATION_HPP

#include "4C_config.hpp"

#include "4C_fracture_coordinate_system.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Fracture
{
  struct NodeContourResult
  {
    /// Configurational-force vector per unit crack-front support.
    Vector3 configurational_force;
    double j_integral;
    Vector3 propagation_direction;
  };

  /**
   * Normalize an integrated configurational force, project scalar J, and determine propagation.
   *
   * The propagation criterion is the normalized configurational force projected onto the plane
   * normal to the crack-front tangent.
   */
  NodeContourResult evaluate_node_contour_result(const Vector3& integrated_configurational_force,
      double front_support_measure, const CoordinateSystem& coordinate_system,
      bool evaluate_propagation_direction = true);
}  // namespace Fracture

FOUR_C_NAMESPACE_CLOSE

#endif
