// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fracture_manager.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fracture_crack_front.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_io_visualization_parameters.hpp"
#include "4C_legacy_enum_definitions_materials.hpp"
#include "4C_linalg_multi_vector.hpp"
#include "4C_solid_ele.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>

FOUR_C_NAMESPACE_OPEN

namespace
{
  using FourC::Fracture::Front;
  using FourC::Fracture::FrontSegment;
  using FourC::Fracture::Manager;
  using FourC::Fracture::Vector3;

  std::vector<int> canonical_line(std::vector<int> line)
  {
    if (line.size() == 3) line = {line[0], line[2], line[1]};
    if (line.front() > line.back()) std::ranges::reverse(line);
    return line;
  }

  std::string sanitized_id(const std::string& id)
  {
    std::string result = id;
    for (char& character : result)
      if (!std::isalnum(static_cast<unsigned char>(character)) && character != '_' &&
          character != '-')
        character = '_';
    return result;
  }

  std::string yaml_quoted(const std::string& value)
  {
    std::string result = "\"";
    for (const char character : value)
    {
      if (character == '"' || character == '\\') result.push_back('\\');
      result.push_back(character);
    }
    result.push_back('"');
    return result;
  }

  void initialize_restarted_ascii_output(
      const std::string& source_file, const std::string& target_file, const double restart_time)
  {
    std::ifstream source(source_file);
    FOUR_C_ASSERT_ALWAYS(
        source.good(), "Cannot read fracture restart output file '{}'.", source_file);

    std::vector<std::string> header;
    std::vector<std::vector<std::string>> steps;
    std::vector<std::string>* current_step = nullptr;
    std::string line;
    while (std::getline(source, line))
    {
      if (line.starts_with("  - step:"))
      {
        steps.emplace_back();
        current_step = &steps.back();
      }
      (current_step == nullptr ? header : *current_step).push_back(line);
    }

    std::ofstream target(target_file, std::ios::trunc);
    FOUR_C_ASSERT_ALWAYS(
        target.good(), "Cannot initialize fracture output file '{}'.", target_file);
    for (const auto& header_line : header) target << header_line << '\n';
    for (const auto& step : steps)
    {
      const auto time_line = std::ranges::find_if(
          step, [](const std::string& value) { return value.starts_with("    time:"); });
      FOUR_C_ASSERT_ALWAYS(time_line != step.end(),
          "Invalid fracture restart output in '{}': step has no time.", source_file);
      const double time = std::stod(time_line->substr(time_line->find(':') + 1));
      if (time > restart_time + 1.0e-14) continue;
      for (const auto& step_line : step) target << step_line << '\n';
    }
  }

  bool contains_quantity(
      const FourC::Fracture::Output& output, const FourC::Fracture::OutputQuantity quantity)
  {
    return std::ranges::find(output.quantities, quantity) != output.quantities.end();
  }

  bool output_contour(const FourC::Fracture::Crack& crack, const int one_based_contour_number)
  {
    if (const auto* preset =
            std::get_if<FourC::Fracture::OutputContourPreset>(&crack.output.contours.value))
      return *preset == FourC::Fracture::OutputContourPreset::all ||
             one_based_contour_number == crack.contours;
    const auto& contours = std::get<std::vector<int>>(crack.output.contours.value);
    return std::ranges::find(contours, one_based_contour_number) != contours.end();
  }

  std::map<int, std::set<int>> node_adjacency(
      const std::set<int>& domain_elements, const std::vector<Manager::ElementTopology>& topology)
  {
    std::map<int, std::set<int>> adjacency;
    for (const auto& element : topology)
    {
      if (!domain_elements.contains(element.volume.element_id)) continue;
      for (std::size_t first = 0; first < element.volume.node_ids.size(); ++first)
      {
        for (std::size_t second = first + 1; second < element.volume.node_ids.size(); ++second)
        {
          adjacency[element.volume.node_ids[first]].insert(element.volume.node_ids[second]);
          adjacency[element.volume.node_ids[second]].insert(element.volume.node_ids[first]);
        }
      }
    }
    return adjacency;
  }

  std::set<int> artificial_boundary_nodes(
      const std::set<int>& domain_elements, const std::vector<Manager::ElementTopology>& topology)
  {
    std::map<std::vector<int>, std::vector<std::pair<int, std::vector<int>>>> faces;
    for (const auto& element : topology)
    {
      for (auto face : element.volume.faces)
      {
        const auto original_face = face;
        std::ranges::sort(face);
        face.erase(std::ranges::unique(face).begin(), face.end());
        if (face.size() < 3) continue;
        faces[face].emplace_back(element.volume.element_id, original_face);
      }
    }

    std::set<int> boundary_nodes;
    for (const auto& [face, incident_elements] : faces)
    {
      if (incident_elements.size() != 2) continue;
      const bool first_inside = domain_elements.contains(incident_elements[0].first);
      const bool second_inside = domain_elements.contains(incident_elements[1].first);
      if (first_inside == second_inside) continue;
      boundary_nodes.insert(face.begin(), face.end());
    }
    return boundary_nodes;
  }

  struct HarmonicExtensionDomain
  {
    std::map<int, std::set<int>> adjacency;
    std::set<int> fixed_nodes;
    std::vector<int> free_nodes;
    std::map<int, std::size_t> free_node_index;
    std::vector<double> diagonal;
    std::vector<std::vector<std::size_t>> free_neighbors;
  };

  HarmonicExtensionDomain make_harmonic_extension_domain(const std::set<int>& front_nodes,
      const std::set<int>& domain_elements, const std::vector<Manager::ElementTopology>& topology)
  {
    HarmonicExtensionDomain domain{.adjacency = node_adjacency(domain_elements, topology),
        .fixed_nodes = artificial_boundary_nodes(domain_elements, topology),
        .free_nodes = {},
        .free_node_index = {},
        .diagonal = {},
        .free_neighbors = {}};
    domain.fixed_nodes.insert(front_nodes.begin(), front_nodes.end());
    for (const auto& [node, neighbors] : domain.adjacency)
    {
      if (domain.fixed_nodes.contains(node)) continue;
      FOUR_C_ASSERT_ALWAYS(
          !neighbors.empty(), "Isolated node {} in fracture contour domain.", node);
      domain.free_node_index.emplace(node, domain.free_nodes.size());
      domain.free_nodes.push_back(node);
    }
    domain.diagonal.resize(domain.free_nodes.size());
    domain.free_neighbors.resize(domain.free_nodes.size());
    for (std::size_t row = 0; row < domain.free_nodes.size(); ++row)
    {
      const auto& neighbors = domain.adjacency.at(domain.free_nodes[row]);
      domain.diagonal[row] = static_cast<double>(neighbors.size());
      for (const int neighbor : neighbors)
        if (const auto index = domain.free_node_index.find(neighbor);
            index != domain.free_node_index.end())
          domain.free_neighbors[row].push_back(index->second);
    }
    return domain;
  }

  std::vector<Vector3> crack_face_interior_vectors(const int front_node, const Front& front,
      const std::vector<Manager::ElementTopology>& topology,
      const std::map<int, Vector3>& coordinates)
  {
    std::set<int> adjacent_front_nodes;
    for (const auto& segment : front.segments)
    {
      const auto position = std::ranges::find(segment.node_ids, front_node);
      if (position == segment.node_ids.end()) continue;
      if (position != segment.node_ids.begin()) adjacent_front_nodes.insert(*(position - 1));
      if (position + 1 != segment.node_ids.end()) adjacent_front_nodes.insert(*(position + 1));
    }

    std::map<std::vector<int>, std::vector<int>> boundary_faces;
    for (const auto& element : topology)
    {
      for (const auto& face : element.volume.faces)
      {
        auto key = face;
        std::ranges::sort(key);
        auto [entry, inserted] = boundary_faces.try_emplace(key, face);
        if (!inserted) entry->second.clear();
      }
    }

    std::vector<Vector3> vectors;
    for (const auto& [key, face] : boundary_faces)
    {
      if (face.empty() || std::ranges::find(face, front_node) == face.end() ||
          !std::ranges::any_of(
              face, [&](const int node) { return adjacent_front_nodes.contains(node); }))
        continue;

      Vector3 vector{};
      int number_of_face_nodes = 0;
      for (const int face_node : face)
      {
        if (face_node == front_node) continue;
        for (int direction = 0; direction < 3; ++direction)
          vector[direction] +=
              coordinates.at(face_node)[direction] - coordinates.at(front_node)[direction];
        ++number_of_face_nodes;
      }
      if (number_of_face_nodes > 0)
      {
        for (double& component : vector) component /= number_of_face_nodes;
        vectors.push_back(vector);
      }
    }
    FOUR_C_ASSERT_ALWAYS(!vectors.empty(),
        "Cannot reconstruct fracture frame at node {}: no crack-face boundary faces were found.",
        front_node);
    return vectors;
  }

  std::map<int, double> harmonic_extension(
      const int target_node, const HarmonicExtensionDomain& domain)
  {
    std::map<int, double> values;
    for (const auto& [node, neighbors] : domain.adjacency)
      values[node] = node == target_node ? 1.0 : 0.0;
    FOUR_C_ASSERT_ALWAYS(values.contains(target_node) && domain.fixed_nodes.contains(target_node),
        "Crack-front node {} is not contained in its contour domain.", target_node);

    if (domain.free_nodes.empty()) return values;

    const auto apply_laplacian = [&](const std::vector<double>& input)
    {
      std::vector<double> output(domain.free_nodes.size(), 0.0);
      for (std::size_t row = 0; row < domain.free_nodes.size(); ++row)
      {
        output[row] = domain.diagonal[row] * input[row];
        for (const std::size_t neighbor : domain.free_neighbors[row])
          output[row] -= input[neighbor];
      }
      return output;
    };
    const auto dot = [](const std::vector<double>& first, const std::vector<double>& second)
    {
      double result = 0.0;
      for (std::size_t i = 0; i < first.size(); ++i) result += first[i] * second[i];
      return result;
    };

    std::vector<double> solution(domain.free_nodes.size(), 0.0);
    std::vector<double> residual(domain.free_nodes.size(), 0.0);
    for (std::size_t row = 0; row < domain.free_nodes.size(); ++row)
      residual[row] = domain.adjacency.at(domain.free_nodes[row]).contains(target_node) ? 1.0 : 0.0;
    std::vector<double> preconditioned_residual(domain.free_nodes.size());
    for (std::size_t i = 0; i < domain.free_nodes.size(); ++i)
      preconditioned_residual[i] = residual[i] / domain.diagonal[i];
    std::vector<double> direction = preconditioned_residual;
    double residual_norm_squared = dot(residual, residual);
    double residual_preconditioned_residual = dot(residual, preconditioned_residual);
    const double tolerance_squared = 1.0e-24 * std::max(1.0, residual_norm_squared);

    const std::size_t maximum_iterations =
        std::max<std::size_t>(100, 10 * domain.free_nodes.size());
    for (std::size_t iteration = 0;
        iteration < maximum_iterations && residual_norm_squared > tolerance_squared; ++iteration)
    {
      const std::vector<double> matrix_direction = apply_laplacian(direction);
      const double denominator = dot(direction, matrix_direction);
      FOUR_C_ASSERT_ALWAYS(
          denominator > 0.0, "The harmonic fracture extension system is not positive definite.");
      const double step = residual_preconditioned_residual / denominator;
      for (std::size_t i = 0; i < domain.free_nodes.size(); ++i)
      {
        solution[i] += step * direction[i];
        residual[i] -= step * matrix_direction[i];
      }
      const double next_residual_norm_squared = dot(residual, residual);
      if (next_residual_norm_squared <= tolerance_squared)
      {
        residual_norm_squared = next_residual_norm_squared;
        break;
      }
      for (std::size_t i = 0; i < domain.free_nodes.size(); ++i)
        preconditioned_residual[i] = residual[i] / domain.diagonal[i];
      const double next_residual_preconditioned_residual = dot(residual, preconditioned_residual);
      const double beta = next_residual_preconditioned_residual / residual_preconditioned_residual;
      for (std::size_t i = 0; i < domain.free_nodes.size(); ++i)
        direction[i] = preconditioned_residual[i] + beta * direction[i];
      residual_norm_squared = next_residual_norm_squared;
      residual_preconditioned_residual = next_residual_preconditioned_residual;
    }
    FOUR_C_ASSERT_ALWAYS(residual_norm_squared <= tolerance_squared,
        "Harmonic fracture extension field did not converge for front node {} after {} iterations.",
        target_node, maximum_iterations);

    for (std::size_t i = 0; i < domain.free_nodes.size(); ++i)
      values.at(domain.free_nodes[i]) = solution[i];
    return values;
  }
}  // namespace

Fracture::Manager::Manager(Core::FE::Discretization& discretization, std::vector<Crack> cracks,
    std::string output_file_base, const Core::IO::VisualizationParameters& visualization_parameters)
    : discretization_(discretization)
{
  topology_ = collect_global_topology();
  coordinates_ = collect_global_coordinates(topology_);

  const bool requests_g = std::ranges::any_of(cracks,
      [](const Crack& crack) { return contains_quantity(crack.output, OutputQuantity::G); });
  if (requests_g)
  {
    int local_unsupported_material = 0;
    for (int local_id = 0; local_id < discretization_.num_my_row_elements(); ++local_id)
    {
      const auto* solid =
          dynamic_cast<const Discret::Elements::Solid<3>*>(discretization_.l_row_element(local_id));
      if (solid == nullptr) continue;
      const auto material_type = solid->solid_material()->material_type();
      if (material_type != Core::Materials::m_stvenant &&
          material_type != Core::Materials::m_orthostvenant)
        local_unsupported_material = 1;
    }
    FOUR_C_ASSERT_ALWAYS(
        Core::Communication::max_all(local_unsupported_material, discretization_.get_comm()) == 0,
        "Fracture output quantity G is supported only for linear-elastic St. Venant-Kirchhoff "
        "materials.");
  }

  std::vector<const Core::Conditions::Condition*> conditions;
  discretization_.get_condition(std::string(crack_front_condition_name), conditions);
  FOUR_C_ASSERT_ALWAYS(conditions.size() == cracks.size(),
      "Found {} CrackFront conditions for {} fracture configurations.", conditions.size(),
      cracks.size());

  std::map<int, const Core::Conditions::Condition*> conditions_by_index;
  for (const auto* condition : conditions)
    conditions_by_index.emplace(condition->parameters().get<int>("configuration_index"), condition);

  fronts_.reserve(cracks.size());
  std::set<std::string> sanitized_crack_ids;
  for (std::size_t i = 0; i < cracks.size(); ++i)
  {
    FOUR_C_ASSERT_ALWAYS(conditions_by_index.contains(static_cast<int>(i)),
        "Missing CrackFront condition for fracture configuration {}.", i);
    fronts_.push_back(build_front_data(
        cracks[i], *conditions_by_index.at(static_cast<int>(i)), topology_, coordinates_, i));
    const std::string output_id = sanitized_id(cracks[i].crack_id);
    FOUR_C_ASSERT_ALWAYS(sanitized_crack_ids.insert(output_id).second,
        "Fracture crack_ids produce duplicate sanitized output name '{}'.", output_id);
    if (cracks[i].output.data_format == DataFormat::ascii)
    {
      fronts_.back().output_file = output_file_base + "_fracture_" + output_id + ".yaml";
      if (Core::Communication::my_mpi_rank(discretization_.get_comm()) == 0)
      {
        if (visualization_parameters.restart_time_ > 0.0)
        {
          const std::string restart_file =
              visualization_parameters.restart_from_name_ + "_fracture_" + output_id + ".yaml";
          initialize_restarted_ascii_output(
              restart_file, fronts_.back().output_file, visualization_parameters.restart_time_);
        }
        else
        {
          std::ofstream output(fronts_.back().output_file, std::ios::trunc);
          FOUR_C_ASSERT_ALWAYS(output.good(), "Cannot initialize fracture output file '{}'.",
              fronts_.back().output_file);
        }
      }
    }
    else
    {
      auto binary_parameters = visualization_parameters;
      binary_parameters.data_format_ = Core::IO::OutputDataFormat::binary;
      fronts_.back().binary_writer = std::make_unique<Core::IO::VisualizationManager>(
          std::move(binary_parameters), discretization_.get_comm(), "fracture_" + output_id);
    }
  }
}

Fracture::Manager::~Manager() = default;

std::vector<Fracture::Manager::ElementTopology> Fracture::Manager::collect_global_topology() const
{
  std::vector<int> serialized;
  for (int local_id = 0; local_id < discretization_.num_my_row_elements(); ++local_id)
  {
    auto* solid =
        dynamic_cast<Discret::Elements::Solid<3>*>(discretization_.l_row_element(local_id));
    if (solid == nullptr) continue;

    const int* node_ids = solid->node_ids();
    const auto face_connectivity = Core::FE::get_ele_node_numbering_surfaces(solid->shape());
    const auto line_connectivity = Core::FE::get_ele_node_numbering_lines(solid->shape());
    serialized.push_back(solid->id());
    serialized.push_back(solid->num_node());
    serialized.insert(serialized.end(), node_ids, node_ids + solid->num_node());
    serialized.push_back(static_cast<int>(face_connectivity.size()));
    for (const auto& face : face_connectivity)
    {
      serialized.push_back(static_cast<int>(face.size()));
      for (const int local_node : face) serialized.push_back(node_ids[local_node]);
    }
    serialized.push_back(static_cast<int>(line_connectivity.size()));
    for (const auto& line : line_connectivity)
    {
      serialized.push_back(static_cast<int>(line.size()));
      for (const int local_node : line) serialized.push_back(node_ids[local_node]);
    }
  }

  const auto gathered = Core::Communication::all_gather(serialized, discretization_.get_comm());
  std::map<int, ElementTopology> topology_by_id;
  for (const auto& rank_data : gathered)
  {
    std::size_t position = 0;
    while (position < rank_data.size())
    {
      ElementTopology element;
      element.volume.element_id = rank_data[position++];
      const int number_of_nodes = rank_data[position++];
      element.volume.node_ids.assign(
          rank_data.begin() + position, rank_data.begin() + position + number_of_nodes);
      position += number_of_nodes;
      const int number_of_faces = rank_data[position++];
      for (int face = 0; face < number_of_faces; ++face)
      {
        const int number_of_face_nodes = rank_data[position++];
        element.volume.faces.emplace_back(
            rank_data.begin() + position, rank_data.begin() + position + number_of_face_nodes);
        position += number_of_face_nodes;
      }
      const int number_of_lines = rank_data[position++];
      for (int line = 0; line < number_of_lines; ++line)
      {
        const int number_of_line_nodes = rank_data[position++];
        std::vector<int> node_ids(
            rank_data.begin() + position, rank_data.begin() + position + number_of_line_nodes);
        position += number_of_line_nodes;
        element.lines.push_back({canonical_line(std::move(node_ids))});
      }
      topology_by_id.emplace(element.volume.element_id, std::move(element));
    }
  }

  std::vector<ElementTopology> topology;
  topology.reserve(topology_by_id.size());
  for (auto& [id, element] : topology_by_id) topology.push_back(std::move(element));
  return topology;
}

std::map<int, Fracture::Vector3> Fracture::Manager::collect_global_coordinates(
    const std::vector<ElementTopology>& topology) const
{
  std::vector<int> local_ids;
  std::vector<double> local_coordinates;
  for (int local_id = 0; local_id < discretization_.num_my_row_nodes(); ++local_id)
  {
    const auto* node = discretization_.l_row_node(local_id);
    local_ids.push_back(node->id());
    local_coordinates.insert(local_coordinates.end(), node->x().begin(), node->x().end());
  }
  const auto gathered_ids = Core::Communication::all_gather(local_ids, discretization_.get_comm());
  const auto gathered_coordinates =
      Core::Communication::all_gather(local_coordinates, discretization_.get_comm());

  std::map<int, Vector3> coordinates;
  for (std::size_t rank = 0; rank < gathered_ids.size(); ++rank)
  {
    FOUR_C_ASSERT_ALWAYS(gathered_coordinates[rank].size() == 3 * gathered_ids[rank].size(),
        "Invalid gathered fracture coordinate data.");
    for (std::size_t node = 0; node < gathered_ids[rank].size(); ++node)
      coordinates[gathered_ids[rank][node]] = {gathered_coordinates[rank][3 * node],
          gathered_coordinates[rank][3 * node + 1], gathered_coordinates[rank][3 * node + 2]};
  }
  return coordinates;
}

Fracture::Manager::FrontData Fracture::Manager::build_front_data(const Crack& crack,
    const Core::Conditions::Condition& condition, const std::vector<ElementTopology>& topology,
    const std::map<int, Vector3>& coordinates, const std::size_t crack_index) const
{
  const std::set<int> front_nodes(condition.get_nodes()->begin(), condition.get_nodes()->end());
  std::set<int> front_neighborhood_nodes;
  for (const auto& element : topology)
    if (std::ranges::any_of(
            element.volume.node_ids, [&](const int node) { return front_nodes.contains(node); }))
      front_neighborhood_nodes.insert(
          element.volume.node_ids.begin(), element.volume.node_ids.end());
  for (const int front_node : front_nodes)
  {
    const auto& front_coordinate = coordinates.at(front_node);
    for (const int other_node : front_neighborhood_nodes)
    {
      if (other_node == front_node) continue;
      const auto& other_coordinate = coordinates.at(other_node);
      double squared_distance = 0.0;
      for (int direction = 0; direction < 3; ++direction)
      {
        const double difference = front_coordinate[direction] - other_coordinate[direction];
        squared_distance += difference * difference;
      }
      FOUR_C_ASSERT_ALWAYS(squared_distance > 1.0e-24,
          "Crack-front node {} of crack_id '{}' is coincident with distinct mesh node {}. "
          "Crack-front nodes must be conforming and shared by all incident solid elements.",
          front_node, crack.crack_id, other_node);
    }
  }

  std::map<std::vector<int>, FrontSegment> front_segments;
  std::vector<VolumeElementTopology> volume_topology;
  volume_topology.reserve(topology.size());
  for (const auto& element : topology)
  {
    volume_topology.push_back(element.volume);
    for (const auto& line : element.lines)
    {
      if (std::ranges::all_of(
              line.node_ids, [&](const int node) { return front_nodes.contains(node); }) &&
          line.node_ids.front() != line.node_ids.back())
        front_segments.emplace(canonical_line(line.node_ids), line);
    }
  }
  FOUR_C_ASSERT_ALWAYS(!front_segments.empty(),
      "No solid-element edges connect the nodes of fracture crack_id '{}'.", crack.crack_id);

  std::vector<FrontSegment> segments;
  for (const auto& [key, segment] : front_segments) segments.push_back(segment);
  Front front = build_front(segments, coordinates);
  std::vector<CoordinateSystem> coordinate_systems;
  coordinate_systems.reserve(front.nodes.size());
  for (const auto& node : front.nodes)
  {
    if (crack.direction || crack.normal)
      coordinate_systems.push_back(
          make_coordinate_system(node.tangent, crack.direction, crack.normal));
    else
      coordinate_systems.push_back(make_coordinate_system_from_geometry(
          node.tangent, crack_face_interior_vectors(node.node_id, front, topology, coordinates)));
  }

  ContourLayers contour_layers = build_contour_layers(volume_topology, front_nodes, crack.contours);
  std::vector<std::vector<std::map<int, double>>> extension_fields(
      front.nodes.size(), std::vector<std::map<int, double>>(crack.contours));
  {
    TEUCHOS_FUNC_TIME_MONITOR("Fracture: construct extension fields");

    for (int contour = 0; contour < crack.contours; ++contour)
    {
      if (contour > 0 && contour_layers.contours[contour] == contour_layers.contours[contour - 1])
      {
        for (std::size_t node = 0; node < front.nodes.size(); ++node)
          extension_fields[node][contour] = extension_fields[node][contour - 1];
        continue;
      }

      const HarmonicExtensionDomain domain =
          make_harmonic_extension_domain(front_nodes, contour_layers.contours[contour], topology);
      std::vector<std::map<int, double>> contour_extension_fields;
      contour_extension_fields.reserve(front.nodes.size());
      for (std::size_t node = 0; node < front.nodes.size(); ++node)
        contour_extension_fields.push_back(harmonic_extension(front.nodes[node].node_id, domain));
      normalize_extension_fields(contour_extension_fields);
      for (std::size_t node = 0; node < front.nodes.size(); ++node)
        extension_fields[node][contour] = std::move(contour_extension_fields[node]);
    }
  }

  return {.crack = crack,
      .front = std::move(front),
      .coordinate_systems = std::move(coordinate_systems),
      .contour_layers = std::move(contour_layers),
      .extension_fields = std::move(extension_fields),
      .results = std::vector<std::vector<NodeContourResult>>(
          front_nodes.size(), std::vector<NodeContourResult>(crack.contours)),
      .output_file = {},
      .binary_writer = nullptr};
}

void Fracture::Manager::evaluate(const double time, const double time_step_size)
{
  TEUCHOS_FUNC_TIME_MONITOR("Fracture: calculate J-integral");

  Teuchos::ParameterList parameters;
  parameters.set("total time", time);
  parameters.set("delta time", time_step_size);

  for (auto& front : fronts_)
  {
    for (std::size_t node = 0; node < front.front.nodes.size(); ++node)
    {
      for (int contour = 0; contour < front.crack.contours; ++contour)
      {
        std::array<double, 3> local_force{0.0, 0.0, 0.0};
        const auto& domain = front.contour_layers.contours[contour];
        const auto& extension = front.extension_fields[node][contour];
        for (int local_id = 0; local_id < discretization_.num_my_row_elements(); ++local_id)
        {
          auto* solid =
              dynamic_cast<Discret::Elements::Solid<3>*>(discretization_.l_row_element(local_id));
          if (solid == nullptr || !domain.contains(solid->id())) continue;

          std::vector<double> nodal_extension(solid->num_node());
          for (int element_node = 0; element_node < solid->num_node(); ++element_node)
            nodal_extension[element_node] = extension.at(solid->node_ids()[element_node]);
          Core::Elements::LocationArray location_array(1);
          solid->location_vector(discretization_, location_array);
          const auto element_force = solid->calculate_configurational_force(
              discretization_, location_array[0].lm_, parameters, nodal_extension);
          for (int direction = 0; direction < 3; ++direction)
            local_force[direction] += element_force[direction];
        }

        Vector3 global_force{};
        for (int direction = 0; direction < 3; ++direction)
          global_force[direction] =
              Core::Communication::sum_all(local_force[direction], discretization_.get_comm());
        front.results[node][contour] = evaluate_node_contour_result(global_force,
            front.front.nodes[node].support_measure, front.coordinate_systems[node],
            contains_quantity(front.crack.output, OutputQuantity::propagation_direction));
      }
    }
  }
}

void Fracture::Manager::write_step(const int step, const double time) const
{
  const int rank = Core::Communication::my_mpi_rank(discretization_.get_comm());
  for (const auto& front : fronts_)
  {
    if (front.crack.output.data_format == DataFormat::binary)
    {
      front.binary_writer->clear_data();
      auto& data = front.binary_writer->get_visualization_data();
      if (rank == 0)
      {
        std::map<int, int> local_point_index;
        for (std::size_t node = 0; node < front.front.nodes.size(); ++node)
        {
          local_point_index.emplace(front.front.nodes[node].node_id, static_cast<int>(node));
          const auto& coordinate = coordinates_.at(front.front.nodes[node].node_id);
          data.get_point_coordinates().insert(
              data.get_point_coordinates().end(), coordinate.begin(), coordinate.end());
        }
        for (const auto& segment : front.front.segments)
        {
          data.get_cell_types().push_back(segment.node_ids.size() == 2 ? 3 : 21);
          data.get_cell_connectivity().push_back(local_point_index.at(segment.node_ids[0]));
          data.get_cell_connectivity().push_back(
              local_point_index.at(segment.node_ids[segment.node_ids.size() == 2 ? 1 : 2]));
          if (segment.node_ids.size() == 3)
            data.get_cell_connectivity().push_back(local_point_index.at(segment.node_ids[1]));
          data.get_cell_offsets().push_back(
              static_cast<Core::IO::index_type>(data.get_cell_connectivity().size()));
        }

        std::vector<int> node_ids;
        std::vector<double> support_measures;
        std::vector<double> directions;
        std::vector<double> normals;
        std::vector<double> tangents;
        for (std::size_t node = 0; node < front.front.nodes.size(); ++node)
        {
          node_ids.push_back(front.front.nodes[node].node_id);
          support_measures.push_back(front.front.nodes[node].support_measure);
          directions.insert(directions.end(),
              front.coordinate_systems[node].propagation_direction.begin(),
              front.coordinate_systems[node].propagation_direction.end());
          normals.insert(normals.end(), front.coordinate_systems[node].normal.begin(),
              front.coordinate_systems[node].normal.end());
          tangents.insert(tangents.end(), front.coordinate_systems[node].tangent.begin(),
              front.coordinate_systems[node].tangent.end());
        }
        data.set_point_data_vector("node_id", node_ids, 1);
        data.set_point_data_vector("support_measure", support_measures, 1);
        data.set_point_data_vector("crack_direction", directions, 3);
        data.set_point_data_vector("crack_normal", normals, 3);
        data.set_point_data_vector("crack_front_tangent", tangents, 3);

        for (int contour = 0; contour < front.crack.contours; ++contour)
        {
          if (!output_contour(front.crack, contour + 1)) continue;
          std::vector<double> j_integrals;
          std::vector<double> propagation_directions;
          for (std::size_t node = 0; node < front.front.nodes.size(); ++node)
          {
            const auto& result = front.results[node][contour];
            j_integrals.push_back(result.j_integral);
            propagation_directions.insert(propagation_directions.end(),
                result.propagation_direction.begin(), result.propagation_direction.end());
          }
          const std::string suffix = "_contour_" + std::to_string(contour + 1);
          if (contains_quantity(front.crack.output, OutputQuantity::J))
            data.set_point_data_vector("J" + suffix, j_integrals, 1);
          if (contains_quantity(front.crack.output, OutputQuantity::G))
            data.set_point_data_vector("G" + suffix, j_integrals, 1);
          if (contains_quantity(front.crack.output, OutputQuantity::propagation_direction))
            data.set_point_data_vector("propagation_direction" + suffix, propagation_directions, 3);
        }
      }
      front.binary_writer->write_to_disk(time, step);
      continue;
    }

    if (rank != 0) continue;
    std::ifstream existing_output(front.output_file);
    const bool new_file =
        !existing_output.good() || existing_output.peek() == std::ifstream::traits_type::eof();
    std::ofstream output(front.output_file, std::ios::app);
    FOUR_C_ASSERT_ALWAYS(
        output.good(), "Cannot open fracture output file '{}'.", front.output_file);
    output << std::setprecision(16);
    if (new_file)
      output << "crack_id: " << yaml_quoted(front.crack.crack_id)
             << "\nsymmetry_scaling: none\nsteps:\n";
    output << "  - step: " << step << "\n    time: " << time << "\n    nodes:\n";
    for (std::size_t node = 0; node < front.front.nodes.size(); ++node)
    {
      output << "      - node_id: " << front.front.nodes[node].node_id << "\n"
             << "        coordinates: [" << coordinates_.at(front.front.nodes[node].node_id)[0]
             << ", " << coordinates_.at(front.front.nodes[node].node_id)[1] << ", "
             << coordinates_.at(front.front.nodes[node].node_id)[2] << "]\n"
             << "        support_measure: " << front.front.nodes[node].support_measure << "\n"
             << "        frame:\n"
             << "          direction: [" << front.coordinate_systems[node].propagation_direction[0]
             << ", " << front.coordinate_systems[node].propagation_direction[1] << ", "
             << front.coordinate_systems[node].propagation_direction[2] << "]\n"
             << "          normal: [" << front.coordinate_systems[node].normal[0] << ", "
             << front.coordinate_systems[node].normal[1] << ", "
             << front.coordinate_systems[node].normal[2] << "]\n"
             << "          tangent: [" << front.coordinate_systems[node].tangent[0] << ", "
             << front.coordinate_systems[node].tangent[1] << ", "
             << front.coordinate_systems[node].tangent[2] << "]\n"
             << "        contours:\n";
      for (int contour = 0; contour < front.crack.contours; ++contour)
      {
        if (!output_contour(front.crack, contour + 1)) continue;
        const auto& result = front.results[node][contour];
        output << "          - contour: " << contour + 1 << "\n";
        if (contains_quantity(front.crack.output, OutputQuantity::J))
          output << "            J: " << result.j_integral << "\n";
        if (contains_quantity(front.crack.output, OutputQuantity::G))
          output << "            G: " << result.j_integral << "\n";
        if (contains_quantity(front.crack.output, OutputQuantity::propagation_direction))
          output << "            propagation_direction: [" << result.propagation_direction[0]
                 << ", " << result.propagation_direction[1] << ", "
                 << result.propagation_direction[2] << "]\n";
      }
    }
  }
}

std::unique_ptr<Core::LinAlg::MultiVector<double>> Fracture::Manager::element_layer_data(
    const std::size_t crack_index) const
{
  FOUR_C_ASSERT_ALWAYS(
      crack_index < fronts_.size(), "Invalid fracture crack index {}.", crack_index);
  auto data = std::make_unique<Core::LinAlg::MultiVector<double>>(
      *discretization_.element_row_map(), 1, true);
  for (const auto& [element_id, layer] : fronts_[crack_index].contour_layers.element_layer)
    if (discretization_.element_row_map()->my_gid(element_id))
      data->replace_global_value(element_id, 0, layer);
  return data;
}

const std::string& Fracture::Manager::crack_id(const std::size_t crack_index) const
{
  FOUR_C_ASSERT_ALWAYS(
      crack_index < fronts_.size(), "Invalid fracture crack index {}.", crack_index);
  return fronts_[crack_index].crack.crack_id;
}

FOUR_C_NAMESPACE_CLOSE
