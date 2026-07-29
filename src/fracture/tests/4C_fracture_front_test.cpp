// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fracture_front.hpp"

#include "4C_fracture_coordinate_system.hpp"
#include "4C_utils_exceptions.hpp"

#include <cmath>

namespace
{
  using namespace FourC;

  void expect_vector_near(
      const Fracture::Vector3& actual, const Fracture::Vector3& expected, const double tolerance)
  {
    for (int i = 0; i < 3; ++i) EXPECT_NEAR(actual[i], expected[i], tolerance);
  }

  TEST(FractureFrontTest, OrdersOpenFrontAndComputesLine2Support)
  {
    const std::map<int, Fracture::Vector3> coordinates = {
        {1, {0.0, 0.0, 0.0}}, {2, {0.0, 0.0, 1.0}}, {3, {0.0, 0.0, 3.0}}};
    const auto front = Fracture::build_front({{{2, 3}}, {{2, 1}}}, coordinates);

    ASSERT_FALSE(front.closed);
    ASSERT_EQ(front.nodes.size(), 3);
    EXPECT_EQ(front.nodes[0].node_id, 1);
    EXPECT_EQ(front.nodes[1].node_id, 2);
    EXPECT_EQ(front.nodes[2].node_id, 3);
    EXPECT_NEAR(front.nodes[0].support_measure, 0.5, 1.0e-14);
    EXPECT_NEAR(front.nodes[1].support_measure, 1.5, 1.0e-14);
    EXPECT_NEAR(front.nodes[2].support_measure, 1.0, 1.0e-14);
    for (const auto& node : front.nodes) expect_vector_near(node.tangent, {0.0, 0.0, 1.0}, 1.0e-14);
  }

  TEST(FractureFrontTest, SupportsQuadraticLineInterpolation)
  {
    const std::map<int, Fracture::Vector3> coordinates = {
        {1, {0.0, 0.0, 0.0}}, {2, {0.0, 0.0, 0.5}}, {3, {0.0, 0.0, 1.0}}};
    const auto front = Fracture::build_front({{{1, 2, 3}}}, coordinates);

    ASSERT_EQ(front.nodes.size(), 3);
    EXPECT_NEAR(front.nodes[0].support_measure, 1.0 / 6.0, 1.0e-14);
    EXPECT_NEAR(front.nodes[1].support_measure, 2.0 / 3.0, 1.0e-14);
    EXPECT_NEAR(front.nodes[2].support_measure, 1.0 / 6.0, 1.0e-14);
  }

  TEST(FractureFrontTest, OrdersClosedFrontDeterministically)
  {
    const std::map<int, Fracture::Vector3> coordinates = {
        {1, {0.0, 0.0, 0.0}}, {2, {1.0, 0.0, 0.0}}, {3, {1.0, 1.0, 0.0}}, {4, {0.0, 1.0, 0.0}}};
    const auto front = Fracture::build_front({{{3, 4}}, {{1, 4}}, {{2, 3}}, {{1, 2}}}, coordinates);

    ASSERT_TRUE(front.closed);
    ASSERT_EQ(front.nodes.size(), 4);
    EXPECT_EQ(front.nodes[0].node_id, 1);
    EXPECT_EQ(front.nodes[1].node_id, 2);
    EXPECT_EQ(front.nodes[2].node_id, 3);
    EXPECT_EQ(front.nodes[3].node_id, 4);
  }

  TEST(FractureFrontTest, RejectsBranchingAndDisconnectedSegments)
  {
    const std::map<int, Fracture::Vector3> coordinates = {{1, {0.0, 0.0, 0.0}},
        {2, {1.0, 0.0, 0.0}}, {3, {2.0, 0.0, 0.0}}, {4, {1.0, 1.0, 0.0}}, {5, {4.0, 0.0, 0.0}},
        {6, {5.0, 0.0, 0.0}}};
    EXPECT_THROW(
        Fracture::build_front({{{1, 2}}, {{2, 3}}, {{2, 4}}}, coordinates), Core::Exception);
    EXPECT_THROW(Fracture::build_front({{{1, 2}}, {{5, 6}}}, coordinates), Core::Exception);
  }

  TEST(FractureCoordinateSystemTest, BuildsFramesFromEitherOrBothInputs)
  {
    const Fracture::Vector3 tangent{0.0, 0.0, 1.0};

    const auto from_direction =
        Fracture::make_coordinate_system(tangent, std::vector<double>{2.0, 0.0, 1.0}, std::nullopt);
    expect_vector_near(from_direction.propagation_direction, {1.0, 0.0, 0.0}, 1.0e-14);
    expect_vector_near(from_direction.normal, {0.0, 1.0, 0.0}, 1.0e-14);

    const auto from_normal =
        Fracture::make_coordinate_system(tangent, std::nullopt, std::vector<double>{0.0, 3.0, 1.0});
    expect_vector_near(from_normal.propagation_direction, {1.0, 0.0, 0.0}, 1.0e-14);
    expect_vector_near(from_normal.normal, {0.0, 1.0, 0.0}, 1.0e-14);

    const auto from_both = Fracture::make_coordinate_system(
        tangent, std::vector<double>{2.0, 0.0, 1.0}, std::vector<double>{0.0, 3.0, 1.0});
    expect_vector_near(from_both.tangent, {0.0, 0.0, 1.0}, 1.0e-14);
  }

  TEST(FractureCoordinateSystemTest, RejectsMissingOrDegenerateInput)
  {
    const Fracture::Vector3 tangent{0.0, 0.0, 1.0};
    EXPECT_THROW(
        Fracture::make_coordinate_system(tangent, std::nullopt, std::nullopt), Core::Exception);
    EXPECT_THROW(
        Fracture::make_coordinate_system(tangent, std::vector<double>{0.0, 0.0, 2.0}, std::nullopt),
        Core::Exception);
  }

  TEST(FractureCoordinateSystemTest, UsesSuppliedFrameToOrientTangent)
  {
    const auto frame = Fracture::make_coordinate_system(
        {0.0, 0.0, 1.0}, std::vector<double>{1.0, 0.0, 0.0}, std::vector<double>{0.0, -1.0, 0.0});
    expect_vector_near(frame.tangent, {0.0, 0.0, -1.0}, 1.0e-14);
  }

  TEST(FractureCoordinateSystemTest, ReconstructsFrameFromCrackFaceGeometry)
  {
    const auto frame = Fracture::make_coordinate_system_from_geometry(
        {0.0, 0.0, 1.0}, {{-2.0, 1.0, 0.5}, {-2.0, -1.0, -0.5}});
    expect_vector_near(frame.propagation_direction, {1.0, 0.0, 0.0}, 1.0e-14);
    expect_vector_near(frame.normal, {0.0, 1.0, 0.0}, 1.0e-14);
    expect_vector_near(frame.tangent, {0.0, 0.0, 1.0}, 1.0e-14);
  }
}  // namespace
