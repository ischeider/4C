// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fracture_input.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

#include <limits>
#include <set>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace
{
  using FourC::Fracture::OutputQuantity;

  std::vector<OutputQuantity> parse_quantities(const std::vector<std::string>& names)
  {
    if (names.empty()) FOUR_C_THROW("Fracture output quantity must not be empty.");

    std::set<std::string> unique_names;
    std::vector<OutputQuantity> quantities;
    quantities.reserve(names.size());

    for (const std::string& name : names)
    {
      if (!unique_names.insert(name).second)
        FOUR_C_THROW("Duplicate fracture output quantity '{}'.", name);

      if (name == "J")
        quantities.push_back(OutputQuantity::J);
      else if (name == "G")
        quantities.push_back(OutputQuantity::G);
      else if (name == "propagation_direction")
        quantities.push_back(OutputQuantity::propagation_direction);
      else
        FOUR_C_THROW("Unsupported fracture output quantity '{}'.", name);
    }

    return quantities;
  }

  void validate_contour_numbers(const std::vector<int>& contour_numbers)
  {
    if (contour_numbers.empty())
      FOUR_C_THROW("Fracture output_contours must not be an empty vector.");

    std::set<int> unique_contours;
    for (const int contour : contour_numbers)
    {
      if (contour <= 0) FOUR_C_THROW("Fracture contour numbers must be positive.");
      if (!unique_contours.insert(contour).second)
        FOUR_C_THROW("Duplicate fracture output contour {}.", contour);
    }
  }

  std::optional<std::vector<double>> get_optional_vector(
      const Teuchos::ParameterList& parameters, const std::string& name)
  {
    if (!parameters.isParameter(name)) return std::nullopt;
    if (parameters.isType<std::optional<std::vector<double>>>(name))
      return parameters.get<std::optional<std::vector<double>>>(name);
    if (parameters.isType<std::vector<double>>(name))
      return parameters.get<std::vector<double>>(name);

    FOUR_C_THROW("Fracture parameter '{}' has an unsupported type.", name);
  }

  double squared_norm(const std::vector<double>& vector)
  {
    return vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
  }

  void validate_vector(const std::optional<std::vector<double>>& vector, const std::string& name)
  {
    if (!vector) return;
    if (vector->size() != 3)
      FOUR_C_THROW("Fracture parameter '{}' must contain exactly three components.", name);
    if (squared_norm(*vector) == 0.0)
      FOUR_C_THROW("Fracture parameter '{}' must not be the zero vector.", name);
  }

  bool are_collinear(const std::vector<double>& first, const std::vector<double>& second)
  {
    const double cross_x = first[1] * second[2] - first[2] * second[1];
    const double cross_y = first[2] * second[0] - first[0] * second[2];
    const double cross_z = first[0] * second[1] - first[1] * second[0];
    const double cross_squared = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
    const double scale = squared_norm(first) * squared_norm(second);
    return cross_squared <= 100.0 * std::numeric_limits<double>::epsilon() * scale;
  }

  FourC::Fracture::OutputContours parse_output_contours(
      const Teuchos::ParameterList& output_parameters, const int number_of_contours)
  {
    if (output_parameters.isType<FourC::Fracture::OutputContourPreset>("output_contours"))
    {
      const auto preset =
          output_parameters.get<FourC::Fracture::OutputContourPreset>("output_contours");
      if (preset != FourC::Fracture::OutputContourPreset::all &&
          preset != FourC::Fracture::OutputContourPreset::outer)
        FOUR_C_THROW("Unsupported fracture output_contours preset.");
      return {.value = preset};
    }

    if (output_parameters.isType<std::vector<int>>("output_contours"))
    {
      const auto& contour_numbers = output_parameters.get<std::vector<int>>("output_contours");
      validate_contour_numbers(contour_numbers);
      for (const int contour : contour_numbers)
      {
        if (contour > number_of_contours)
          FOUR_C_THROW("Fracture output contour {} exceeds the configured {} contours.", contour,
              number_of_contours);
      }
      return {.value = contour_numbers};
    }

    FOUR_C_THROW(
        "Fracture output_contours must be 'all', 'outer', or a vector of contour numbers.");
  }
}  // namespace

Core::IO::InputSpec Fracture::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;

  return list("STRUCTURAL DYNAMIC/FRACTURE",
      all_of({
          parameter<std::string>("crack_id",
              {.description = "External mesh node-set name identifying the crack front."}),
          parameter<std::optional<std::vector<double>>>("direction",
              {.description = "Optional crack propagation direction; normalization is deferred to "
                              "geometry processing.",
                  .size = 3}),
          parameter<std::optional<std::vector<double>>>("normal",
              {.description = "Optional crack-plane normal; normalization is deferred to geometry "
                              "processing.",
                  .size = 3}),
          parameter<int>("contours", {.description = "Number of positive contour layers.",
                                         .validator = Validators::positive<int>()}),
          parameter<PropagationCriterion>("propagation_criterion",
              {.description = "Crack propagation based on maximum configurational force, not MERR.",
                  .default_value = PropagationCriterion::maximum_configurational_force}),
          group("output",
              {
                  parameter<std::vector<std::string>>("quantity",
                      {.description =
                              "Nonempty unique selection of J, G, and propagation_direction. J is "
                              "not automatically doubled for symmetry.",
                          .on_parse_callback =
                              [](Core::IO::InputParameterContainer& container)
                          {
                            parse_quantities(container.get<std::vector<std::string>>("quantity"));
                          }}),
                  one_of({
                      parameter<OutputContourPreset>("output_contours",
                          {.description = "Output all contours or only the outer contour."}),
                      parameter<std::vector<int>>("output_contours",
                          {.description = "Unique positive contour numbers to output.",
                              .on_parse_callback =
                                  [](Core::IO::InputParameterContainer& container)
                              {
                                validate_contour_numbers(
                                    container.get<std::vector<int>>("output_contours"));
                              }}),
                  }),
                  parameter<DataFormat>(
                      "data_format", {.description = "Fracture output data format.",
                                         .default_value = DataFormat::ascii}),
              }),
      }),
      {.description = "Optional fracture mechanics definitions for structural dynamics.",
          .required = false});
}

std::vector<Fracture::Crack> Fracture::parse_parameters(
    const std::vector<Teuchos::ParameterList>& parameters)
{
  std::set<std::string> crack_ids;
  std::vector<Crack> cracks;
  cracks.reserve(parameters.size());

  for (const Teuchos::ParameterList& parameters_for_crack : parameters)
  {
    Crack crack;
    crack.crack_id = parameters_for_crack.get<std::string>("crack_id");
    if (crack.crack_id.empty()) FOUR_C_THROW("Fracture crack_id must not be empty.");
    if (!crack_ids.insert(crack.crack_id).second)
      FOUR_C_THROW("Duplicate fracture crack_id '{}'.", crack.crack_id);

    crack.direction = get_optional_vector(parameters_for_crack, "direction");
    crack.normal = get_optional_vector(parameters_for_crack, "normal");
    validate_vector(crack.direction, "direction");
    validate_vector(crack.normal, "normal");
    if (crack.direction && crack.normal && are_collinear(*crack.direction, *crack.normal))
      FOUR_C_THROW("Fracture direction and normal must not be collinear.");

    crack.contours = parameters_for_crack.get<int>("contours");
    if (crack.contours <= 0) FOUR_C_THROW("Fracture contours must be positive.");

    crack.propagation_criterion =
        parameters_for_crack.get<PropagationCriterion>("propagation_criterion");
    if (crack.propagation_criterion != PropagationCriterion::maximum_configurational_force)
      FOUR_C_THROW("Unsupported fracture propagation criterion.");

    const Teuchos::ParameterList& output_parameters = parameters_for_crack.sublist("output");
    crack.output.quantities =
        parse_quantities(output_parameters.get<std::vector<std::string>>("quantity"));
    crack.output.contours = parse_output_contours(output_parameters, crack.contours);
    crack.output.data_format = output_parameters.get<DataFormat>("data_format");
    if (crack.output.data_format != DataFormat::ascii &&
        crack.output.data_format != DataFormat::binary)
      FOUR_C_THROW("Unsupported fracture output data format.");

    cracks.push_back(std::move(crack));
  }

  return cracks;
}

FOUR_C_NAMESPACE_CLOSE
