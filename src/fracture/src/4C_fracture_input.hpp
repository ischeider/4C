// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FRACTURE_INPUT_HPP
#define FOUR_C_FRACTURE_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Fracture
{
  enum class PropagationCriterion
  {
    maximum_configurational_force
  };

  enum class OutputQuantity
  {
    J,
    G,
    propagation_direction
  };

  enum class OutputContourPreset
  {
    all,
    outer
  };

  enum class DataFormat
  {
    ascii,
    binary
  };

  struct OutputContours
  {
    std::variant<OutputContourPreset, std::vector<int>> value;
  };

  struct Output
  {
    std::vector<OutputQuantity> quantities;
    OutputContours contours;
    DataFormat data_format = DataFormat::ascii;
  };

  struct Crack
  {
    std::string crack_id;
    std::optional<std::vector<double>> direction;
    std::optional<std::vector<double>> normal;
    int contours;
    PropagationCriterion propagation_criterion =
        PropagationCriterion::maximum_configurational_force;
    Output output;
  };

  /**
   * Input specification for the optional STRUCTURAL DYNAMIC/FRACTURE list.
   *
   * The propagation criterion is the maximum configurational force. It is not the maximum energy
   * release rate (MERR). Reported J values are not automatically doubled for symmetry.
   */
  Core::IO::InputSpec valid_parameters();

  /// Convert the ParameterList representation produced by valid_parameters() to validated data.
  std::vector<Crack> parse_parameters(const std::vector<Teuchos::ParameterList>& parameters);
}  // namespace Fracture

FOUR_C_NAMESPACE_CLOSE

#endif
