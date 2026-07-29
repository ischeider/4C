// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FRACTURE_MANAGER_HPP
#define FOUR_C_FRACTURE_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_fracture_contours.hpp"
#include "4C_fracture_coordinate_system.hpp"
#include "4C_fracture_front.hpp"
#include "4C_fracture_input.hpp"
#include "4C_fracture_propagation.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}

namespace Core::Conditions
{
  class Condition;
}

namespace Core::IO
{
  struct VisualizationParameters;
  class VisualizationManager;
}  // namespace Core::IO

namespace Core::LinAlg
{
  template <typename T>
  class MultiVector;
}

namespace Fracture
{
  class Manager
  {
   public:
    struct ElementTopology
    {
      VolumeElementTopology volume;
      std::vector<FrontSegment> lines;
    };

    Manager(Core::FE::Discretization& discretization, std::vector<Crack> cracks,
        std::string output_file_base,
        const Core::IO::VisualizationParameters& visualization_parameters);
    ~Manager();

    void evaluate(double time, double time_step_size);

    void write_step(int step, double time) const;

    [[nodiscard]] std::unique_ptr<Core::LinAlg::MultiVector<double>> element_layer_data(
        std::size_t crack_index) const;

    [[nodiscard]] std::size_t number_of_cracks() const { return fronts_.size(); }

    [[nodiscard]] const std::string& crack_id(std::size_t crack_index) const;

   private:
    struct FrontData
    {
      Crack crack;
      Front front;
      std::vector<CoordinateSystem> coordinate_systems;
      ContourLayers contour_layers;
      std::vector<std::vector<std::map<int, double>>> extension_fields;
      std::vector<std::vector<NodeContourResult>> results;
      std::string output_file;
      std::unique_ptr<Core::IO::VisualizationManager> binary_writer;
    };

    std::vector<ElementTopology> collect_global_topology() const;

    std::map<int, Vector3> collect_global_coordinates(
        const std::vector<ElementTopology>& topology) const;

    FrontData build_front_data(const Crack& crack, const Core::Conditions::Condition& condition,
        const std::vector<ElementTopology>& topology, const std::map<int, Vector3>& coordinates,
        std::size_t crack_index) const;

    Core::FE::Discretization& discretization_;
    std::vector<ElementTopology> topology_;
    std::map<int, Vector3> coordinates_;
    std::vector<FrontData> fronts_;
  };
}  // namespace Fracture

FOUR_C_NAMESPACE_CLOSE

#endif
