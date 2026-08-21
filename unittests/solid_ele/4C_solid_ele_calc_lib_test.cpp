// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_solid_ele_calc_lib.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_linalg_tensor_symmetric_einstein.hpp"
#include "4C_solid_ele_calc_eas_helpers.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <numbers>
#include <type_traits>

namespace
{
  using namespace FourC;

  TEST(EvaluateParameterCoordinateCentroid, DisTypeHex)
  {
    // only tested for hex8, but equivalent for hex18, hex27, ...
    const auto distype = Core::FE::CellType::hex8;

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> xi_centroid_ref{};

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> xi_centroid =
        Discret::Elements::evaluate_parameter_coordinate_centroid<distype>();

    FOUR_C_EXPECT_NEAR(xi_centroid, xi_centroid_ref, 1e-15);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypeTet)
  {
    // only tested for tet4, but equivalent for tet10
    const auto distype = Core::FE::CellType::tet4;

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> xi_centroid_ref = {
        {0.25, 0.25, 0.25}};

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> xi_centroid =
        Discret::Elements::evaluate_parameter_coordinate_centroid<distype>();

    FOUR_C_EXPECT_NEAR(xi_centroid, xi_centroid_ref, 1e-15);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypePyramid)
  {
    const auto distype = Core::FE::CellType::pyramid5;

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> xi_centroid_ref = {
        {0.0, 0.0, 0.25}};

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> xi_centroid =
        Discret::Elements::evaluate_parameter_coordinate_centroid<distype>();

    FOUR_C_EXPECT_NEAR(xi_centroid, xi_centroid_ref, 1e-15);
  }

  TEST(EvaluateParameterCoordinateCentroid, DisTypeWedge)
  {
    const auto distype = Core::FE::CellType::wedge6;

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> xi_centroid_ref = {
        {1.0 / 3.0, 1.0 / 3.0, 0.0}};


    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> xi_centroid =
        Discret::Elements::evaluate_parameter_coordinate_centroid<distype>();

    FOUR_C_EXPECT_NEAR(xi_centroid, xi_centroid_ref, 1e-15);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeHex)
  {
    const auto distype = Core::FE::CellType::hex8;

    Core::LinAlg::Matrix<Discret::Elements::Internal::num_dim<distype>, 1>
        reference_coords_centroid_ref(Core::LinAlg::Initialization::zero);

    Discret::Elements::ElementNodes<distype> nodal_coordinates;

    nodal_coordinates.reference_coordinates(0, 0) = 0;
    nodal_coordinates.reference_coordinates(1, 0) = 0;
    nodal_coordinates.reference_coordinates(2, 0) = 0;

    nodal_coordinates.reference_coordinates(0, 1) = 4;
    nodal_coordinates.reference_coordinates(1, 1) = 0;
    nodal_coordinates.reference_coordinates(2, 1) = 0;

    nodal_coordinates.reference_coordinates(0, 2) = 4;
    nodal_coordinates.reference_coordinates(1, 2) = 1;
    nodal_coordinates.reference_coordinates(2, 2) = 0;

    nodal_coordinates.reference_coordinates(0, 3) = 0;
    nodal_coordinates.reference_coordinates(1, 3) = 1;
    nodal_coordinates.reference_coordinates(2, 3) = 0;

    nodal_coordinates.reference_coordinates(0, 4) = 0;
    nodal_coordinates.reference_coordinates(1, 4) = 0;
    nodal_coordinates.reference_coordinates(2, 4) = 2;

    nodal_coordinates.reference_coordinates(0, 5) = 4;
    nodal_coordinates.reference_coordinates(1, 5) = 0;
    nodal_coordinates.reference_coordinates(2, 5) = 2;

    nodal_coordinates.reference_coordinates(0, 6) = 4;
    nodal_coordinates.reference_coordinates(1, 6) = 1;
    nodal_coordinates.reference_coordinates(2, 6) = 2;

    nodal_coordinates.reference_coordinates(0, 7) = 0;
    nodal_coordinates.reference_coordinates(1, 7) = 1;
    nodal_coordinates.reference_coordinates(2, 7) = 2;


    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> x_centroid_ref = {
        {2, 0.5, 1.0}};

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> x_centroid =
        Discret::Elements::evaluate_reference_coordinate_centroid<distype>(nodal_coordinates);

    FOUR_C_EXPECT_NEAR(x_centroid, x_centroid_ref, 1e-15);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeTet)
  {
    const auto distype = Core::FE::CellType::tet4;

    Core::LinAlg::Matrix<Discret::Elements::Internal::num_dim<distype>, 1>
        reference_coords_centroid_ref(Core::LinAlg::Initialization::zero);

    Discret::Elements::ElementNodes<distype> nodal_coordinates;

    nodal_coordinates.reference_coordinates(0, 0) = 0;
    nodal_coordinates.reference_coordinates(1, 0) = 0;
    nodal_coordinates.reference_coordinates(2, 0) = 0;

    nodal_coordinates.reference_coordinates(0, 1) = 1;
    nodal_coordinates.reference_coordinates(1, 1) = 0;
    nodal_coordinates.reference_coordinates(2, 1) = 0;

    nodal_coordinates.reference_coordinates(0, 2) = 0;
    nodal_coordinates.reference_coordinates(1, 2) = 2;
    nodal_coordinates.reference_coordinates(2, 2) = 0;

    nodal_coordinates.reference_coordinates(0, 3) = 0;
    nodal_coordinates.reference_coordinates(1, 3) = 0;
    nodal_coordinates.reference_coordinates(2, 3) = 4;

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> x_centroid_ref = {
        {0.25, 0.5, 1.0}};

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> x_centroid =
        Discret::Elements::evaluate_reference_coordinate_centroid<distype>(nodal_coordinates);

    FOUR_C_EXPECT_NEAR(x_centroid, x_centroid_ref, 1e-15);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypePyramid)
  {
    const auto distype = Core::FE::CellType::pyramid5;

    Core::LinAlg::Matrix<Discret::Elements::Internal::num_dim<distype>, 1>
        reference_coords_centroid_ref(Core::LinAlg::Initialization::zero);

    Discret::Elements::ElementNodes<distype> nodal_coordinates;

    nodal_coordinates.reference_coordinates(0, 0) = -2;
    nodal_coordinates.reference_coordinates(1, 0) = -1;
    nodal_coordinates.reference_coordinates(2, 0) = 0;

    nodal_coordinates.reference_coordinates(0, 1) = 4;
    nodal_coordinates.reference_coordinates(1, 1) = -1;
    nodal_coordinates.reference_coordinates(2, 1) = 0;

    nodal_coordinates.reference_coordinates(0, 2) = 4;
    nodal_coordinates.reference_coordinates(1, 2) = 1;
    nodal_coordinates.reference_coordinates(2, 2) = 0;

    nodal_coordinates.reference_coordinates(0, 3) = -2;
    nodal_coordinates.reference_coordinates(1, 3) = 1;
    nodal_coordinates.reference_coordinates(2, 3) = 0;

    nodal_coordinates.reference_coordinates(0, 4) = 1;
    nodal_coordinates.reference_coordinates(1, 4) = 0;
    nodal_coordinates.reference_coordinates(2, 4) = 4;


    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> x_centroid_ref = {
        {1.0, 0.0, 1.0}};

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> x_centroid =
        Discret::Elements::evaluate_reference_coordinate_centroid<distype>(nodal_coordinates);

    FOUR_C_EXPECT_NEAR(x_centroid, x_centroid_ref, 1e-15);
  }

  TEST(EvaluateReferenceCoordinateCentroid, DisTypeWedge)
  {
    const auto distype = Core::FE::CellType::wedge6;

    Core::LinAlg::Matrix<Discret::Elements::Internal::num_dim<distype>, 1>
        reference_coords_centroid_ref(Core::LinAlg::Initialization::zero);

    Discret::Elements::ElementNodes<distype> nodal_coordinates;

    nodal_coordinates.reference_coordinates(0, 0) = 0;
    nodal_coordinates.reference_coordinates(1, 0) = 0;
    nodal_coordinates.reference_coordinates(2, 0) = 0;

    nodal_coordinates.reference_coordinates(0, 1) = 3;
    nodal_coordinates.reference_coordinates(1, 1) = 0;
    nodal_coordinates.reference_coordinates(2, 1) = 0;

    nodal_coordinates.reference_coordinates(0, 2) = 0;
    nodal_coordinates.reference_coordinates(1, 2) = 6;
    nodal_coordinates.reference_coordinates(2, 2) = 0;

    nodal_coordinates.reference_coordinates(0, 3) = 0;
    nodal_coordinates.reference_coordinates(1, 3) = 0;
    nodal_coordinates.reference_coordinates(2, 3) = 1;

    nodal_coordinates.reference_coordinates(0, 4) = 3;
    nodal_coordinates.reference_coordinates(1, 4) = 0;
    nodal_coordinates.reference_coordinates(2, 4) = 1;

    nodal_coordinates.reference_coordinates(0, 5) = 0;
    nodal_coordinates.reference_coordinates(1, 5) = 6;
    nodal_coordinates.reference_coordinates(2, 5) = 1;

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> x_centroid_ref = {
        {1.0, 2.0, 0.5}};

    Core::LinAlg::Tensor<double, Discret::Elements::Internal::num_dim<distype>> x_centroid =
        Discret::Elements::evaluate_reference_coordinate_centroid<distype>(nodal_coordinates);

    FOUR_C_EXPECT_NEAR(x_centroid, x_centroid_ref, 1e-15);
  }

  namespace
  {
    template <Core::FE::CellType celltype>
    Discret::Elements::ElementNodes<celltype> get_element_nodes()
    {
      std::srand(0);
      Discret::Elements::ElementNodes<celltype> nodes;
      Core::LinAlg::SerialDenseMatrix reference_nodes =
          Core::FE::get_ele_node_numbering_nodes_paramspace(celltype);

      for (int i = 0; i < Core::FE::num_nodes(celltype); ++i)
      {
        for (int d = 0; d < Core::FE::dim<celltype>; ++d)
        {
          nodes.reference_coordinates(d, i) =
              reference_nodes(d, i) + 0.03 * ((double)std::rand()) / RAND_MAX;
          nodes.displacements(d, i) = 0.02 * ((double)std::rand()) / RAND_MAX;


          nodes.current_coordinates(d, i) =
              nodes.reference_coordinates(d, i) + nodes.displacements(d, i);
        }
      }

      return nodes;
    }

    template <Core::FE::CellType celltype>
    Discret::Elements::JacobianMapping<celltype> get_jacobian_mapping(
        const Discret::Elements::ElementNodes<celltype>& nodes)
    {
      return Discret::Elements::evaluate_jacobian_mapping(
          Discret::Elements::evaluate_shape_functions_and_derivs({{0.1, 0.2, 0.3}}, nodes), nodes);
    }

    template <Core::FE::CellType celltype>
    Discret::Elements::Stress<celltype> get_stress()
    {
      Discret::Elements::Stress<celltype> stress;
      stress.pk2_ = Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{{
          {1.1, 1.2, 1.3},
          {1.2, 2.2, 2.3},
          {1.3, 2.3, 3.3},
      }});
      stress.cmat_ = Core::LinAlg::einsum_sym<"AC", "BD">(stress.pk2_, stress.pk2_) +
                     Core::LinAlg::einsum_sym<"AD", "BC">(stress.pk2_, stress.pk2_);
      return stress;
    }

    template <Core::FE::CellType... celltypes>
    using make_celltyped_test =
        testing::Types<std::integral_constant<Core::FE::CellType, celltypes>...>;
  }  // namespace

  using celltype_list = make_celltyped_test<Core::FE::CellType::hex8, Core::FE::CellType::hex18,
      Core::FE::CellType::hex20, Core::FE::CellType::hex27, Core::FE::CellType::tet4,
      Core::FE::CellType::tet10, Core::FE::CellType::pyramid5, Core::FE::CellType::wedge6,
      Core::FE::CellType::wedge15>;

  template <typename>
  struct SolidEleCalcLibTest : public testing::Test
  {
  };
  TYPED_TEST_SUITE(SolidEleCalcLibTest, celltype_list);

  TYPED_TEST(SolidEleCalcLibTest, BFreeElementStiffnessMatrix)
  {
    constexpr Core::FE::CellType celltype = TypeParam();
    constexpr int num_nodes = Core::FE::num_nodes(celltype);
    constexpr int num_dim = Core::FE::dim<celltype>;
    constexpr int num_dofs = num_nodes * num_dim;
    constexpr int num_str = num_dim * (num_dim + 1) / 2;

    const Discret::Elements::ElementNodes<celltype> nodes = get_element_nodes<celltype>();
    const FourC::Discret::Elements::JacobianMapping<celltype> jacobian_mapping =
        get_jacobian_mapping(nodes);
    const Discret::Elements::SpatialMaterialMapping<celltype> spatial_material_mapping =
        Discret::Elements::evaluate_spatial_material_mapping(jacobian_mapping, nodes);
    Discret::Elements::Stress<celltype> stress = get_stress<celltype>();
    const double integration_factor = 1.0;


    // compute stiffness matrix using the B-operator
    Core::LinAlg::Matrix<num_dofs, num_dofs> stiffness_matrix_bop{};
    Core::LinAlg::Matrix<num_str, num_dofs> Bop =
        Discret::Elements::evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping);

    add_elastic_stiffness_matrix(Bop, stress, integration_factor, stiffness_matrix_bop);
    add_geometric_stiffness_matrix(
        jacobian_mapping, stress.pk2_, integration_factor, stiffness_matrix_bop);


    // compute stiffness matrix without the B-operator
    Core::LinAlg::Matrix<num_dofs, num_dofs> stiffness_matrix_bfree{};
    Discret::Elements::add_stiffness_matrix(jacobian_mapping,
        spatial_material_mapping.deformation_gradient_, stress, integration_factor,
        stiffness_matrix_bfree);

    // Compare both for equality
    FOUR_C_EXPECT_NEAR(stiffness_matrix_bop, stiffness_matrix_bfree, 1e-13);
  }

  namespace
  {
    template <Core::FE::CellType celltype>
    Discret::Elements::ElementNodes<celltype> make_axisymmetric_nodes()
    {
      Discret::Elements::ElementNodes<celltype> nodes;
      const Core::LinAlg::SerialDenseMatrix parameter_space_nodes =
          Core::FE::get_ele_node_numbering_nodes_paramspace(celltype);

      for (int i = 0; i < Core::FE::num_nodes(celltype); ++i)
      {
        nodes.reference_coordinates(0, i) = parameter_space_nodes(0, i) + 2.0;
        nodes.reference_coordinates(1, i) = parameter_space_nodes(1, i) + 1.0;
        nodes.displacements(0, i) = 0.08 * nodes.reference_coordinates(0, i);
        nodes.displacements(1, i) = -0.03 * nodes.reference_coordinates(1, i);
      }
      nodes.current_coordinates = nodes.reference_coordinates;
      nodes.current_coordinates.update(1.0, nodes.displacements, 1.0);
      return nodes;
    }

    Discret::Elements::ElementNodes<Core::FE::CellType::quad4> make_axisymmetric_quad4_nodes(
        const double inner_radius)
    {
      constexpr auto celltype = Core::FE::CellType::quad4;
      Discret::Elements::ElementNodes<celltype> nodes;
      nodes.reference_coordinates(0, 0) = inner_radius;
      nodes.reference_coordinates(0, 1) = inner_radius + 1.0;
      nodes.reference_coordinates(0, 2) = inner_radius + 1.0;
      nodes.reference_coordinates(0, 3) = inner_radius;
      nodes.reference_coordinates(1, 0) = 0.0;
      nodes.reference_coordinates(1, 1) = 0.0;
      nodes.reference_coordinates(1, 2) = 1.0;
      nodes.reference_coordinates(1, 3) = 1.0;

      constexpr double radial_stretch = 0.08;
      constexpr double axial_stretch = -0.03;
      for (int i = 0; i < 4; ++i)
      {
        nodes.displacements(0, i) = radial_stretch * nodes.reference_coordinates(0, i);
        nodes.displacements(1, i) = axial_stretch * nodes.reference_coordinates(1, i);
      }
      nodes.current_coordinates = nodes.reference_coordinates;
      nodes.current_coordinates.update(1.0, nodes.displacements, 1.0);
      return nodes;
    }

    Discret::Elements::Stress<Core::FE::CellType::quad4> make_stvenant_kirchhoff_stress(
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& strain)
    {
      constexpr double lambda = 2.3;
      constexpr double mu = 1.7;
      Discret::Elements::Stress<Core::FE::CellType::quad4> stress;
      stress.pk2_3d_ = lambda * Core::LinAlg::trace(strain) *
                           Core::LinAlg::TensorGenerators::identity<double, 3, 3> +
                       2.0 * mu * strain;
      stress.cmat_3d_ =
          2.0 * mu * Core::LinAlg::TensorGenerators::symmetric_identity<double, 3, 3, 3, 3> +
          lambda * Core::LinAlg::dyadic(Core::LinAlg::TensorGenerators::identity<double, 3, 3>,
                       Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
      return stress;
    }
  }  // namespace

  using axisymmetric_celltype_list =
      make_celltyped_test<Core::FE::CellType::quad4, Core::FE::CellType::quad8,
          Core::FE::CellType::quad9, Core::FE::CellType::tri3, Core::FE::CellType::tri6>;

  template <typename>
  struct AxisymmetricSolidCellTypeTest : public testing::Test
  {
  };
  TYPED_TEST_SUITE(AxisymmetricSolidCellTypeTest, axisymmetric_celltype_list);

  TYPED_TEST(AxisymmetricSolidCellTypeTest, KinematicsAndThicknessIndependentWeight)
  {
    constexpr Core::FE::CellType celltype = TypeParam();
    const auto nodes = make_axisymmetric_nodes<celltype>();
    const auto xi = Discret::Elements::evaluate_parameter_coordinate_centroid<celltype>();
    const auto shape_functions =
        Discret::Elements::evaluate_shape_functions_and_derivs<celltype>(xi, nodes);
    const auto jacobian_mapping =
        Discret::Elements::evaluate_jacobian_mapping(shape_functions, nodes);
    const auto kinematics = Discret::Elements::evaluate_axisymmetric_kinematics(
        nodes, shape_functions, jacobian_mapping, false);

    EXPECT_NEAR(kinematics.deformation_gradient(0, 0), 1.08, 1e-13);
    EXPECT_NEAR(kinematics.deformation_gradient(1, 1), 0.97, 1e-13);
    EXPECT_NEAR(kinematics.deformation_gradient(2, 2), 1.08, 1e-13);
    EXPECT_NEAR(kinematics.deformation_gradient(0, 1), 0.0, 1e-13);
    EXPECT_NEAR(kinematics.deformation_gradient(1, 0), 0.0, 1e-13);

    const auto integration = Core::FE::create_gauss_integration<celltype>(
        Discret::Elements::get_gauss_rule_stiffness_matrix<celltype>());
    auto integrated_volume = [&](const double thickness)
    {
      const Discret::Elements::ElementProperties<celltype> properties{
          .reference_thickness = thickness,
          .plane_assumption = Discret::Elements::PlaneAssumption::axisymmetric};
      double volume = 0.0;
      Discret::Elements::for_each_gauss_point(nodes, properties, integration,
          [&](const auto& xi, const auto& shape_functions, const auto& jacobian_mapping,
              const double integration_factor, const int gp) { volume += integration_factor; });
      return volume;
    };

    EXPECT_GT(integrated_volume(1.0), 0.0);
    EXPECT_NEAR(integrated_volume(1.0), integrated_volume(99.0), 1e-13);
  }

  TEST(AxisymmetricSolid, Kinematics)
  {
    constexpr auto celltype = Core::FE::CellType::quad4;
    const auto nodes = make_axisymmetric_quad4_nodes(1.0);
    const Core::LinAlg::Tensor<double, 2> xi{};
    const auto shape_functions =
        Discret::Elements::evaluate_shape_functions_and_derivs<celltype>(xi, nodes);
    const auto jacobian_mapping =
        Discret::Elements::evaluate_jacobian_mapping(shape_functions, nodes);

    const auto linear = Discret::Elements::evaluate_axisymmetric_kinematics(
        nodes, shape_functions, jacobian_mapping, true);
    EXPECT_NEAR(linear.green_lagrange_strain(0, 0), 0.08, 1e-14);
    EXPECT_NEAR(linear.green_lagrange_strain(1, 1), -0.03, 1e-14);
    EXPECT_NEAR(linear.green_lagrange_strain(2, 2), 0.08, 1e-14);
    EXPECT_NEAR(linear.green_lagrange_strain(0, 1), 0.0, 1e-14);

    const auto nonlinear = Discret::Elements::evaluate_axisymmetric_kinematics(
        nodes, shape_functions, jacobian_mapping, false);
    EXPECT_NEAR(nonlinear.deformation_gradient(0, 0), 1.08, 1e-14);
    EXPECT_NEAR(nonlinear.deformation_gradient(1, 1), 0.97, 1e-14);
    EXPECT_NEAR(nonlinear.deformation_gradient(2, 2), 1.08, 1e-14);
    EXPECT_NEAR(nonlinear.green_lagrange_strain(0, 0), 0.5 * (1.08 * 1.08 - 1.0), 1e-14);
    EXPECT_NEAR(nonlinear.green_lagrange_strain(1, 1), 0.5 * (0.97 * 0.97 - 1.0), 1e-14);
    EXPECT_NEAR(nonlinear.green_lagrange_strain(2, 2), 0.5 * (1.08 * 1.08 - 1.0), 1e-14);
  }

  TEST(AxisymmetricSolid, EasEnhancedStrainDerivatives)
  {
    Core::LinAlg::Matrix<3, 4> m_tilde(Core::LinAlg::Initialization::zero);
    m_tilde(0, 0) = 1.2;
    m_tilde(1, 1) = 2.3;
    m_tilde(2, 2) = 3.4;
    m_tilde(2, 3) = 4.5;

    std::array<Core::LinAlg::SymmetricTensor<double, 3, 3>, 4> enhanced_strain_derivatives;
    for (int parameter = 0; parameter < 4; ++parameter)
      enhanced_strain_derivatives[parameter] =
          Discret::Elements::evaluate_axisymmetric_eas_strain_derivative<Core::FE::CellType::quad4,
              Discret::Elements::EasType::eastype_q4_4>(m_tilde, parameter);

    EXPECT_DOUBLE_EQ(enhanced_strain_derivatives[0](0, 0), 1.2);
    EXPECT_DOUBLE_EQ(enhanced_strain_derivatives[1](1, 1), 2.3);
    EXPECT_DOUBLE_EQ(enhanced_strain_derivatives[2](0, 1), 1.7);
    EXPECT_DOUBLE_EQ(enhanced_strain_derivatives[3](0, 1), 2.25);
    for (const auto& derivative : enhanced_strain_derivatives)
    {
      EXPECT_DOUBLE_EQ(derivative(2, 2), 0.0);
      EXPECT_DOUBLE_EQ(derivative(0, 2), 0.0);
      EXPECT_DOUBLE_EQ(derivative(1, 2), 0.0);
    }
  }

  TEST(AxisymmetricSolid, EasKinematicsPreservesHoopStrain)
  {
    constexpr auto celltype = Core::FE::CellType::quad4;
    constexpr auto eastype = Discret::Elements::EasType::eastype_q4_4;
    const auto nodes = make_axisymmetric_quad4_nodes(1.0);
    const Core::LinAlg::Tensor<double, 2> xi{{0.25, -0.3}};
    const auto shape_functions =
        Discret::Elements::evaluate_shape_functions_and_derivs<celltype>(xi, nodes);
    const auto jacobian_mapping =
        Discret::Elements::evaluate_jacobian_mapping(shape_functions, nodes);
    const auto centroid_transformation =
        Discret::Elements::evaluate_centroid_transformation<celltype>(nodes);
    Discret::Elements::EasIterationData<celltype, eastype> eas_iteration_data;
    eas_iteration_data.alpha(0) = 0.04;
    eas_iteration_data.alpha(1) = -0.02;
    eas_iteration_data.alpha(2) = 0.03;
    eas_iteration_data.alpha(3) = -0.01;

    const auto kinematics =
        Discret::Elements::evaluate_axisymmetric_eas_kinematics<celltype, eastype>(nodes,
            centroid_transformation, xi, shape_functions, jacobian_mapping, eas_iteration_data);

    EXPECT_DOUBLE_EQ(
        kinematics.enhanced_gl(2, 2), kinematics.displacement_based.green_lagrange_strain(2, 2));
    EXPECT_NE(
        kinematics.enhanced_gl(0, 0), kinematics.displacement_based.green_lagrange_strain(0, 0));
    EXPECT_NE(
        kinematics.enhanced_gl(1, 1), kinematics.displacement_based.green_lagrange_strain(1, 1));
    EXPECT_NE(
        kinematics.enhanced_gl(0, 1), kinematics.displacement_based.green_lagrange_strain(0, 1));
  }

  TEST(AxisymmetricSolid, EasCondensationIncludesHoopCoupling)
  {
    constexpr auto celltype = Core::FE::CellType::quad4;
    constexpr auto eastype = Discret::Elements::EasType::eastype_q4_4;
    const auto nodes = make_axisymmetric_quad4_nodes(1.0);
    const Core::LinAlg::Tensor<double, 2> xi{{0.25, -0.3}};
    const auto shape_functions =
        Discret::Elements::evaluate_shape_functions_and_derivs<celltype>(xi, nodes);
    const auto jacobian_mapping =
        Discret::Elements::evaluate_jacobian_mapping(shape_functions, nodes);
    const auto centroid_transformation =
        Discret::Elements::evaluate_centroid_transformation<celltype>(nodes);
    Discret::Elements::EasIterationData<celltype, eastype> eas_iteration_data;
    const auto kinematics =
        Discret::Elements::evaluate_axisymmetric_eas_kinematics<celltype, eastype>(nodes,
            centroid_transformation, xi, shape_functions, jacobian_mapping, eas_iteration_data);
    const auto stress = make_stvenant_kirchhoff_stress(kinematics.enhanced_gl);
    auto stress_without_hoop_coupling = stress;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l)
            if (i == 2 || j == 2 || k == 2 || l == 2)
              stress_without_hoop_coupling.cmat_3d_(i, j, k, l) = 0.0;

    Discret::Elements::EasIterationData<celltype, eastype> data_with_hoop_coupling;
    Discret::Elements::EasIterationData<celltype, eastype> data_without_hoop_coupling;
    Discret::Elements::integrate_axisymmetric_eas(
        kinematics, shape_functions, jacobian_mapping, stress, 1.0, data_with_hoop_coupling);
    Discret::Elements::integrate_axisymmetric_eas(kinematics, shape_functions, jacobian_mapping,
        stress_without_hoop_coupling, 1.0, data_without_hoop_coupling);

    bool radial_coupling_found = false;
    for (int node = 0; node < 4; ++node)
      for (int parameter = 0; parameter < 4; ++parameter)
      {
        radial_coupling_found |=
            std::abs(data_with_hoop_coupling.Kda(2 * node, parameter) -
                     data_without_hoop_coupling.Kda(2 * node, parameter)) > 1.0e-12;
        EXPECT_NEAR(data_with_hoop_coupling.Kda(2 * node + 1, parameter),
            data_without_hoop_coupling.Kda(2 * node + 1, parameter), 1.0e-14);
      }
    EXPECT_TRUE(radial_coupling_found);
  }

  TEST(AxisymmetricSolid, AxisLimit)
  {
    constexpr auto celltype = Core::FE::CellType::quad4;
    const auto nodes = make_axisymmetric_quad4_nodes(0.0);
    const Core::LinAlg::Tensor<double, 2> xi{{-1.0, 0.0}};
    const auto shape_functions =
        Discret::Elements::evaluate_shape_functions_and_derivs<celltype>(xi, nodes);
    const auto jacobian_mapping =
        Discret::Elements::evaluate_jacobian_mapping(shape_functions, nodes);

    const auto kinematics = Discret::Elements::evaluate_axisymmetric_kinematics(
        nodes, shape_functions, jacobian_mapping, false);
    EXPECT_TRUE(kinematics.use_axis_limit);
    EXPECT_NEAR(kinematics.deformation_gradient(2, 2), 1.08, 1e-14);
  }

  TEST(AxisymmetricSolid, ReferenceVolumeWeight)
  {
    constexpr auto celltype = Core::FE::CellType::quad4;
    const auto nodes = make_axisymmetric_quad4_nodes(1.0);
    const auto integration = Core::FE::create_gauss_integration<celltype>(
        Discret::Elements::get_gauss_rule_stiffness_matrix<celltype>());
    const Discret::Elements::ElementProperties<celltype> properties{.reference_thickness = 99.0,
        .plane_assumption = Discret::Elements::PlaneAssumption::axisymmetric};

    double volume = 0.0;
    Discret::Elements::for_each_gauss_point(nodes, properties, integration,
        [&](const auto& xi, const auto& shape_functions, const auto& jacobian_mapping,
            const double integration_factor, const int gp) { volume += integration_factor; });
    EXPECT_NEAR(volume, 3.0 * std::numbers::pi, 1e-13);
  }

  TEST(AxisymmetricSolid, ConsistentFiniteStrainTangent)
  {
    constexpr auto celltype = Core::FE::CellType::quad4;
    constexpr int num_dofs = 8;
    auto nodes = make_axisymmetric_quad4_nodes(1.0);
    nodes.displacements(0, 2) += 0.02;
    nodes.displacements(1, 3) -= 0.01;
    nodes.current_coordinates = nodes.reference_coordinates;
    nodes.current_coordinates.update(1.0, nodes.displacements, 1.0);
    const Core::LinAlg::Tensor<double, 2> xi{{0.13, -0.21}};

    auto evaluate_force = [&](const auto& evaluated_nodes)
    {
      const auto shape_functions =
          Discret::Elements::evaluate_shape_functions_and_derivs<celltype>(xi, evaluated_nodes);
      const auto jacobian_mapping =
          Discret::Elements::evaluate_jacobian_mapping(shape_functions, evaluated_nodes);
      const auto kinematics = Discret::Elements::evaluate_axisymmetric_kinematics(
          evaluated_nodes, shape_functions, jacobian_mapping, false);
      const auto stress = make_stvenant_kirchhoff_stress(kinematics.green_lagrange_strain);
      Core::LinAlg::Matrix<num_dofs, 1> force{};
      Discret::Elements::add_axisymmetric_internal_force_vector(
          kinematics, shape_functions, jacobian_mapping, stress, 1.0, false, force);
      return force;
    };

    const auto shape_functions =
        Discret::Elements::evaluate_shape_functions_and_derivs<celltype>(xi, nodes);
    const auto jacobian_mapping =
        Discret::Elements::evaluate_jacobian_mapping(shape_functions, nodes);
    const auto kinematics = Discret::Elements::evaluate_axisymmetric_kinematics(
        nodes, shape_functions, jacobian_mapping, false);
    const auto stress = make_stvenant_kirchhoff_stress(kinematics.green_lagrange_strain);
    Core::LinAlg::Matrix<num_dofs, num_dofs> stiffness{};
    Discret::Elements::add_axisymmetric_stiffness_matrix(
        kinematics, shape_functions, jacobian_mapping, stress, 1.0, false, stiffness);

    constexpr double perturbation = 1e-7;
    for (int column = 0; column < num_dofs; ++column)
    {
      auto nodes_plus = nodes;
      auto nodes_minus = nodes;
      const int node = column / 2;
      const int component = column % 2;
      nodes_plus.displacements(component, node) += perturbation;
      nodes_minus.displacements(component, node) -= perturbation;
      nodes_plus.current_coordinates = nodes_plus.reference_coordinates;
      nodes_plus.current_coordinates.update(1.0, nodes_plus.displacements, 1.0);
      nodes_minus.current_coordinates = nodes_minus.reference_coordinates;
      nodes_minus.current_coordinates.update(1.0, nodes_minus.displacements, 1.0);
      const auto force_plus = evaluate_force(nodes_plus);
      const auto force_minus = evaluate_force(nodes_minus);
      for (int row = 0; row < num_dofs; ++row)
      {
        const double finite_difference =
            (force_plus(row) - force_minus(row)) / (2.0 * perturbation);
        EXPECT_NEAR(stiffness(row, column), finite_difference, 2e-8);
      }
    }
  }
}  // namespace