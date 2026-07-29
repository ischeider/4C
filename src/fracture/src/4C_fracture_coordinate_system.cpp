// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fracture_coordinate_system.hpp"

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

  Vector3 cross(const Vector3& left, const Vector3& right)
  {
    return {left[1] * right[2] - left[2] * right[1], left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0]};
  }

  Vector3 scale(const Vector3& vector, const double factor)
  {
    return {factor * vector[0], factor * vector[1], factor * vector[2]};
  }

  Vector3 subtract(const Vector3& left, const Vector3& right)
  {
    return {left[0] - right[0], left[1] - right[1], left[2] - right[2]};
  }

  Vector3 normalize(const Vector3& vector, const std::string_view name)
  {
    const double length = std::sqrt(dot(vector, vector));
    FOUR_C_ASSERT_ALWAYS(length > 100.0 * std::numeric_limits<double>::epsilon(),
        "Cannot construct fracture coordinate system: {} is parallel to the crack-front tangent.",
        name);
    return scale(vector, 1.0 / length);
  }

  Vector3 from_input(const std::vector<double>& vector)
  {
    FOUR_C_ASSERT_ALWAYS(vector.size() == 3, "Fracture coordinate vectors must have size three.");
    return {vector[0], vector[1], vector[2]};
  }

  Vector3 project_normal_to(const Vector3& vector, const Vector3& axis)
  {
    return subtract(vector, scale(axis, dot(vector, axis)));
  }
}  // namespace

Fracture::CoordinateSystem Fracture::make_coordinate_system(const Vector3& front_tangent,
    const std::optional<std::vector<double>>& direction,
    const std::optional<std::vector<double>>& normal)
{
  FOUR_C_ASSERT_ALWAYS(direction || normal,
      "Cannot reconstruct fracture direction and normal without crack-face geometry. Provide at "
      "least one of direction or normal.");

  Vector3 z = normalize(front_tangent, "front tangent");
  if (direction && normal)
  {
    const Vector3 x = normalize(project_normal_to(from_input(*direction), z), "direction");
    const Vector3 y =
        normalize(project_normal_to(project_normal_to(from_input(*normal), z), x), "normal");
    const Vector3 oriented_z = cross(x, y);
    return {.propagation_direction = x, .normal = y, .tangent = oriented_z};
  }

  if (direction)
  {
    const Vector3 x = normalize(project_normal_to(from_input(*direction), z), "direction");
    return {.propagation_direction = x, .normal = cross(z, x), .tangent = z};
  }

  const Vector3 y = normalize(project_normal_to(from_input(*normal), z), "normal");
  return {.propagation_direction = cross(y, z), .normal = y, .tangent = z};
}

Fracture::CoordinateSystem Fracture::make_coordinate_system_from_geometry(
    const Vector3& front_tangent, const std::vector<Vector3>& crack_face_interior_vectors)
{
  const Vector3 z = normalize(front_tangent, "front tangent");
  Vector3 crack_face_interior{};
  for (const auto& vector : crack_face_interior_vectors)
  {
    const Vector3 projected = project_normal_to(vector, z);
    for (int direction = 0; direction < 3; ++direction)
      crack_face_interior[direction] += projected[direction];
  }
  const Vector3 x = scale(normalize(crack_face_interior, "crack-face geometry"), -1.0);
  return {.propagation_direction = x, .normal = cross(z, x), .tangent = z};
}

FOUR_C_NAMESPACE_CLOSE
