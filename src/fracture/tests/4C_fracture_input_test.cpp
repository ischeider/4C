// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fracture_input.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_io_yaml.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

namespace
{
  using namespace FourC;

  Teuchos::ParameterList match(const std::string& yaml)
  {
    Core::IO::InputParameterContainer container;
    ryml::Tree tree = Core::IO::init_yaml_tree_with_exceptions();
    ryml::NodeRef root = tree.rootref();
    ryml::parse_in_arena(ryml::to_csubstr(yaml), root);
    Fracture::valid_parameters().match(Core::IO::ConstYamlNodeRef(root, ""), container);

    Teuchos::ParameterList parameters;
    container.to_teuchos_parameter_list(parameters);
    return parameters;
  }

  std::vector<Fracture::Crack> match_and_parse(const std::string& yaml)
  {
    const Teuchos::ParameterList parameters = match(yaml);
    return Fracture::parse_parameters(
        parameters.get<std::vector<Teuchos::ParameterList>>("STRUCTURAL DYNAMIC/FRACTURE"));
  }

  std::string match_and_emit(const std::string& yaml)
  {
    Core::IO::InputParameterContainer container;
    ryml::Tree input_tree = Core::IO::init_yaml_tree_with_exceptions();
    ryml::NodeRef input_root = input_tree.rootref();
    ryml::parse_in_arena(ryml::to_csubstr(yaml), input_root);

    const auto spec = Fracture::valid_parameters();
    spec.match(Core::IO::ConstYamlNodeRef(input_root, ""), container);

    ryml::Tree output_tree;
    output_tree.rootref() |= ryml::MAP;
    spec.emit(Core::IO::YamlNodeRef(output_tree.rootref(), ""), container);
    return ryml::emitrs_yaml<std::string>(output_tree);
  }

  TEST(FractureInputTest, ParsesValidInputAndDefaults)
  {
    const auto cracks = match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 3
    output:
      quantity: [J, G, propagation_direction]
      output_contours: all
  - crack_id: second_front
    direction: [2.0, 0.0, 0.0]
    normal: [0.0, 3.0, 0.0]
    contours: 4
    propagation_criterion: maximum_configurational_force
    output:
      quantity: [J]
      output_contours: [1, 4]
      data_format: binary
)");

    ASSERT_EQ(cracks.size(), 2);
    EXPECT_EQ(cracks[0].crack_id, "front");
    EXPECT_EQ(cracks[0].direction, std::nullopt);
    EXPECT_EQ(cracks[0].normal, std::nullopt);
    EXPECT_EQ(cracks[0].propagation_criterion,
        Fracture::PropagationCriterion::maximum_configurational_force);
    EXPECT_EQ(cracks[0].output.data_format, Fracture::DataFormat::ascii);
    EXPECT_EQ(std::get<Fracture::OutputContourPreset>(cracks[0].output.contours.value),
        Fracture::OutputContourPreset::all);

    EXPECT_EQ(cracks[1].direction, (std::vector<double>{2.0, 0.0, 0.0}));
    EXPECT_EQ(cracks[1].normal, (std::vector<double>{0.0, 3.0, 0.0}));
    EXPECT_EQ(cracks[1].output.data_format, Fracture::DataFormat::binary);
    EXPECT_EQ(
        std::get<std::vector<int>>(cracks[1].output.contours.value), (std::vector<int>{1, 4}));
  }

  TEST(FractureInputTest, RejectsInvalidVectors)
  {
    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    direction: [0.0, 0.0, 0.0]
    contours: 1
    output: {quantity: [J], output_contours: outer}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    direction: [1.0, 0.0, 0.0]
    normal: [-2.0, 0.0, 0.0]
    contours: 1
    output: {quantity: [J], output_contours: outer}
)"),
        Core::Exception);
  }

  TEST(FractureInputTest, RejectsDuplicateCrackIds)
  {
    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 1
    output: {quantity: [J], output_contours: all}
  - crack_id: front
    contours: 2
    output: {quantity: [G], output_contours: outer}
)"),
        Core::Exception);
  }

  TEST(FractureInputTest, RejectsInvalidContourSelections)
  {
    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 2
    output: {quantity: [J], output_contours: [1, 3]}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 2
    output: {quantity: [J], output_contours: [1, 1]}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 2
    output: {quantity: [J], output_contours: [0, 1]}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 0
    output: {quantity: [J], output_contours: all}
)"),
        Core::Exception);
  }

  TEST(FractureInputTest, RejectsInvalidQuantitiesAndSelections)
  {
    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 1
    output: {quantity: [J, J], output_contours: all}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 1
    output: {quantity: [], output_contours: all}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 1
    output: {quantity: [MERR], output_contours: all}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 1
    propagation_criterion: MERR
    output: {quantity: [J], output_contours: all}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 1
    output: {quantity: [J], output_contours: middle}
)"),
        Core::Exception);

    EXPECT_THROW(match_and_parse(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 1
    output: {quantity: [J], output_contours: all, data_format: csv}
)"),
        Core::Exception);
  }

  TEST(FractureInputTest, EmitsBothOutputContourRepresentations)
  {
    const auto presets = match_and_parse(match_and_emit(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 2
    output: {quantity: [J], output_contours: outer}
)"));
    ASSERT_EQ(presets.size(), 1);
    EXPECT_EQ(std::get<Fracture::OutputContourPreset>(presets[0].output.contours.value),
        Fracture::OutputContourPreset::outer);

    const auto numbers = match_and_parse(match_and_emit(R"(
STRUCTURAL DYNAMIC/FRACTURE:
  - crack_id: front
    contours: 2
    output: {quantity: [J], output_contours: [1, 2]}
)"));
    ASSERT_EQ(numbers.size(), 1);
    EXPECT_EQ(
        std::get<std::vector<int>>(numbers[0].output.contours.value), (std::vector<int>{1, 2}));
  }
}  // namespace
