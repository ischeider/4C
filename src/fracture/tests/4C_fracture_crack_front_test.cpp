// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fracture_crack_front.hpp"

#include "4C_fem_condition.hpp"
#include "4C_fracture_input.hpp"
#include "4C_utils_exceptions.hpp"

namespace
{
  using namespace FourC;

  Fracture::Crack crack(std::string crack_id)
  {
    return Fracture::Crack{.crack_id = std::move(crack_id), .contours = 1};
  }

  TEST(FractureCrackFrontTest, ResolvesNamedNodeSet)
  {
    Core::Conditions::InputNodeSets node_sets;
    node_sets.node_sets = {{17, {8, 3, 5}}};
    node_sets.node_sets_names = {{"front", {17}}};

    auto conditions = Fracture::make_crack_front_conditions({crack("front")}, node_sets);

    ASSERT_EQ(conditions.size(), 1);
    EXPECT_EQ(conditions[0]->node_set_name(), "front");
    EXPECT_EQ(*conditions[0]->get_nodes(), (std::vector<int>{3, 5, 8}));
    EXPECT_EQ(conditions[0]->parameters().get<int>("configuration_index"), 0);
    EXPECT_EQ(conditions[0]->parameters().get<std::string>("crack_id"), "front");
  }

  TEST(FractureCrackFrontTest, ResolvesLegacyDline)
  {
    Core::Conditions::InputNodeSets node_sets;
    node_sets.dline_fenode = {{1, 4}, {6, 2}};

    auto conditions = Fracture::make_crack_front_conditions({crack("2")}, node_sets);

    ASSERT_EQ(conditions.size(), 1);
    EXPECT_EQ(conditions[0]->id(), 1);
    EXPECT_EQ(*conditions[0]->get_nodes(), (std::vector<int>{2, 6}));
  }

  TEST(FractureCrackFrontTest, NamedNodeSetTakesPrecedenceOverNumericLegacyId)
  {
    Core::Conditions::InputNodeSets node_sets;
    node_sets.dline_fenode = {{1, 2}};
    node_sets.node_sets = {{9, {7, 8}}};
    node_sets.node_sets_names = {{"1", {9}}};

    auto conditions = Fracture::make_crack_front_conditions({crack("1")}, node_sets);

    ASSERT_EQ(conditions.size(), 1);
    EXPECT_EQ(conditions[0]->node_set_name(), "1");
    EXPECT_EQ(*conditions[0]->get_nodes(), (std::vector<int>{7, 8}));
  }

  TEST(FractureCrackFrontTest, RejectsMissingAndAmbiguousNodeSets)
  {
    Core::Conditions::InputNodeSets node_sets;
    EXPECT_THROW(
        Fracture::make_crack_front_conditions({crack("missing")}, node_sets), Core::Exception);

    node_sets.node_sets = {{4, {1, 2}}, {5, {3, 4}}};
    node_sets.node_sets_names = {{"front", {4, 5}}};
    EXPECT_THROW(
        Fracture::make_crack_front_conditions({crack("front")}, node_sets), Core::Exception);
  }
}  // namespace
