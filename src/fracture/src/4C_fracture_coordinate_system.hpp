// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FRACTURE_COORDINATE_SYSTEM_HPP
#define FOUR_C_FRACTURE_COORDINATE_SYSTEM_HPP

#include "4C_config.hpp"

#include "4C_fracture_front.hpp"

#include <optional>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Fracture
{
  struct CoordinateSystem
  {
    Vector3 propagation_direction;
    Vector3 normal;
    Vector3 tangent;
  };

  /**
   * Construct x=propagation, y=normal, z=tangent as a right-handed orthonormal frame.
   *
   * At least one user direction is required until crack-face geometry reconstruction is available.
   */
  CoordinateSystem make_coordinate_system(const Vector3& front_tangent,
      const std::optional<std::vector<double>>& direction,
      const std::optional<std::vector<double>>& normal);

  /**
   * Construct a local frame when neither direction nor normal is prescribed.
   *
   * The supplied vectors point from the crack-front node into adjacent crack-face boundary faces.
   * Their component normal to the front tangent points opposite to crack advance.
   */
  CoordinateSystem make_coordinate_system_from_geometry(
      const Vector3& front_tangent, const std::vector<Vector3>& crack_face_interior_vectors);
}  // namespace Fracture

FOUR_C_NAMESPACE_CLOSE

#endif
