// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fracture_propagation.hpp"

#include "4C_utils_exceptions.hpp"

#include <cmath>
#include <limits>

FOUR_C_NAMESPACE_OPEN

namespace
{
  using FourC::Fracture::Vector3;

  double dot(const Vector3& left, const Vector3& right)
  {
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
  }

  Vector3 scale(const Vector3& vector, const double factor)
  {
    return {factor * vector[0], factor * vector[1], factor * vector[2]};
  }

  Vector3 subtract(const Vector3& left, const Vector3& right)
  {
    return {left[0] - right[0], left[1] - right[1], left[2] - right[2]};
  }
}  // namespace

Fracture::NodeContourResult Fracture::evaluate_node_contour_result(
    const Vector3& integrated_configurational_force, const double front_support_measure,
    const CoordinateSystem& coordinate_system, const bool evaluate_propagation_direction)
{
  FOUR_C_ASSERT_ALWAYS(front_support_measure > 0.0,
      "Crack-front support measure must be positive for configurational-force normalization.");

  const Vector3 force = scale(integrated_configurational_force, 1.0 / front_support_measure);
  const double j_integral = dot(force, coordinate_system.propagation_direction);
  const Vector3 in_plane_force =
      subtract(force, scale(coordinate_system.tangent, dot(force, coordinate_system.tangent)));
  const double in_plane_norm = std::sqrt(dot(in_plane_force, in_plane_force));
  if (!evaluate_propagation_direction)
    return {.configurational_force = force, .j_integral = j_integral, .propagation_direction = {}};

  FOUR_C_ASSERT_ALWAYS(in_plane_norm > 100.0 * std::numeric_limits<double>::epsilon(),
      "Maximum configurational-force direction is undefined because the in-plane force vanishes.");
  return {.configurational_force = force,
      .j_integral = j_integral,
      .propagation_direction = scale(in_plane_force, 1.0 / in_plane_norm)};
}

FOUR_C_NAMESPACE_CLOSE
