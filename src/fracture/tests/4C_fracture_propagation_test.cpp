// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fracture_propagation.hpp"

#include "4C_utils_exceptions.hpp"

#include <cmath>

namespace
{
  using namespace FourC;

  TEST(FracturePropagationTest, ProjectsJAndMaximumConfigurationalForceDirection)
  {
    const Fracture::CoordinateSystem frame{.propagation_direction = {1.0, 0.0, 0.0},
        .normal = {0.0, 1.0, 0.0},
        .tangent = {0.0, 0.0, 1.0}};

    const auto result = Fracture::evaluate_node_contour_result({6.0, 8.0, 20.0}, 2.0, frame);

    EXPECT_EQ(result.configurational_force, (Fracture::Vector3{3.0, 4.0, 10.0}));
    EXPECT_DOUBLE_EQ(result.j_integral, 3.0);
    EXPECT_NEAR(result.propagation_direction[0], 3.0 / 5.0, 1.0e-14);
    EXPECT_NEAR(result.propagation_direction[1], 4.0 / 5.0, 1.0e-14);
    EXPECT_NEAR(result.propagation_direction[2], 0.0, 1.0e-14);
  }

  TEST(FracturePropagationTest, SupportsBackwardAndRejectsUndefinedDirection)
  {
    const Fracture::CoordinateSystem frame{.propagation_direction = {1.0, 0.0, 0.0},
        .normal = {0.0, 1.0, 0.0},
        .tangent = {0.0, 0.0, 1.0}};

    EXPECT_THROW(
        Fracture::evaluate_node_contour_result({0.0, 0.0, 1.0}, 1.0, frame), Core::Exception);
    const auto backward = Fracture::evaluate_node_contour_result({-1.0, 1.0, 0.0}, 1.0, frame);
    EXPECT_DOUBLE_EQ(backward.j_integral, -1.0);
    EXPECT_NEAR(backward.propagation_direction[0], -std::sqrt(0.5), 1.0e-14);
    EXPECT_NEAR(backward.propagation_direction[1], std::sqrt(0.5), 1.0e-14);
    EXPECT_THROW(
        Fracture::evaluate_node_contour_result({1.0, 0.0, 0.0}, 0.0, frame), Core::Exception);
  }

  TEST(FracturePropagationTest, AllowsZeroForceWhenDirectionIsNotRequested)
  {
    const Fracture::CoordinateSystem frame{.propagation_direction = {1.0, 0.0, 0.0},
        .normal = {0.0, 1.0, 0.0},
        .tangent = {0.0, 0.0, 1.0}};
    const auto result = Fracture::evaluate_node_contour_result({0.0, 0.0, 0.0}, 1.0, frame, false);
    EXPECT_DOUBLE_EQ(result.j_integral, 0.0);
    EXPECT_EQ(result.propagation_direction, (Fracture::Vector3{}));
  }
}  // namespace
