// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fracture_contours.hpp"

#include "4C_utils_exceptions.hpp"

namespace
{
  using namespace FourC;

  Fracture::VolumeElementTopology hex(const int id, const std::array<int, 8>& nodes)
  {
    return {.element_id = id,
        .node_ids = {nodes.begin(), nodes.end()},
        .faces = {{nodes[0], nodes[1], nodes[2], nodes[3]},
            {nodes[4], nodes[5], nodes[6], nodes[7]}, {nodes[0], nodes[1], nodes[5], nodes[4]},
            {nodes[1], nodes[2], nodes[6], nodes[5]}, {nodes[2], nodes[3], nodes[7], nodes[6]},
            {nodes[3], nodes[0], nodes[4], nodes[7]}}};
  }

  TEST(FractureContoursTest, BuildsNodeAdjacentLayers)
  {
    const std::vector<Fracture::VolumeElementTopology> elements = {
        hex(10, {1, 2, 3, 4, 5, 6, 7, 8}), hex(20, {2, 9, 10, 3, 6, 11, 12, 7}),
        hex(30, {9, 13, 14, 10, 11, 15, 16, 12}),
        // Shares only node 1 with element 10.
        hex(40, {1, 17, 18, 19, 20, 21, 22, 23}), hex(50, {24, 25, 26, 27, 28, 29, 30, 31})};

    const auto layers = Fracture::build_contour_layers(elements, {4}, 3);

    EXPECT_EQ(layers.element_layer.at(10), 1);
    EXPECT_EQ(layers.element_layer.at(20), 2);
    EXPECT_EQ(layers.element_layer.at(30), 3);
    EXPECT_EQ(layers.element_layer.at(40), 2);
    EXPECT_EQ(layers.element_layer.at(50), 0);
    EXPECT_EQ(layers.contours[0], (std::set<int>{10}));
    EXPECT_EQ(layers.contours[1], (std::set<int>{10, 20, 40}));
    EXPECT_EQ(layers.contours[2], (std::set<int>{10, 20, 30, 40}));
  }

  TEST(FractureContoursTest, SupportsTruncatedOpenDomains)
  {
    const std::vector<Fracture::VolumeElementTopology> elements = {
        hex(1, {1, 2, 3, 4, 5, 6, 7, 8}), hex(2, {2, 9, 10, 3, 6, 11, 12, 7})};

    const auto layers = Fracture::build_contour_layers(elements, {1, 4}, 4);

    ASSERT_EQ(layers.contours.size(), 4);
    EXPECT_EQ(layers.contours[0], (std::set<int>{1}));
    EXPECT_EQ(layers.contours[1], (std::set<int>{1, 2}));
    EXPECT_EQ(layers.contours[2], layers.contours[1]);
    EXPECT_EQ(layers.contours[3], layers.contours[1]);
  }

  TEST(FractureContoursTest, SupportsCollapsedCrackTipElements)
  {
    const std::vector<Fracture::VolumeElementTopology> elements = {
        hex(1, {1, 2, 3, 1, 5, 6, 7, 5}), hex(2, {2, 8, 9, 3, 6, 10, 11, 7})};

    const auto layers = Fracture::build_contour_layers(elements, {1, 5}, 2);

    EXPECT_EQ(layers.element_layer.at(1), 1);
    EXPECT_EQ(layers.element_layer.at(2), 2);
    EXPECT_EQ(layers.contours[1], (std::set<int>{1, 2}));
  }

  TEST(FractureContoursTest, NormalizesExtensionFieldsToPartitionOfUnity)
  {
    std::vector<std::map<int, double>> fields = {
        {{1, 1.0}, {2, 0.0}, {3, 0.2}, {4, 0.0}}, {{1, 0.0}, {2, 1.0}, {3, 0.3}, {4, 0.0}}};

    Fracture::normalize_extension_fields(fields);

    EXPECT_DOUBLE_EQ(fields[0].at(1), 1.0);
    EXPECT_DOUBLE_EQ(fields[1].at(2), 1.0);
    EXPECT_DOUBLE_EQ(fields[0].at(3), 0.4);
    EXPECT_DOUBLE_EQ(fields[1].at(3), 0.6);
    EXPECT_DOUBLE_EQ(fields[0].at(4), 0.0);
    EXPECT_DOUBLE_EQ(fields[1].at(4), 0.0);
  }

  TEST(FractureContoursTest, NormalizesSingleExtensionFieldToInteriorPlateau)
  {
    std::vector<std::map<int, double>> fields = {{{1, 1.0}, {2, 0.25}, {3, 0.0}}};

    Fracture::normalize_extension_fields(fields);

    EXPECT_DOUBLE_EQ(fields[0].at(1), 1.0);
    EXPECT_DOUBLE_EQ(fields[0].at(2), 1.0);
    EXPECT_DOUBLE_EQ(fields[0].at(3), 0.0);
  }

  TEST(FractureContoursTest, RejectsIncompatibleExtensionFields)
  {
    std::vector<std::map<int, double>> empty_fields;
    EXPECT_THROW(Fracture::normalize_extension_fields(empty_fields), Core::Exception);

    std::vector<std::map<int, double>> incompatible_fields = {
        {{1, 1.0}, {2, 0.0}}, {{1, 1.0}, {3, 0.0}}};
    EXPECT_THROW(Fracture::normalize_extension_fields(incompatible_fields), Core::Exception);
  }

  TEST(FractureContoursTest, RejectsInvalidTopology)
  {
    EXPECT_THROW(Fracture::build_contour_layers({}, {1}, 1), Core::Exception);

    auto first = hex(1, {1, 2, 3, 4, 5, 6, 7, 8});
    auto second = first;
    second.element_id = 2;
    auto third = first;
    third.element_id = 3;
    EXPECT_THROW(Fracture::build_contour_layers({first, second, third}, {1}, 1), Core::Exception);
  }
}  // namespace
