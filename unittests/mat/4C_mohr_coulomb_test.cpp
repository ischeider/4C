// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_comm_pack_buffer.hpp"
#include "4C_global_data.hpp"
#include "4C_io_input_parameter_container.templates.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_plasticmohrcoulomb.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <array>
#include <cmath>
#include <numbers>

namespace
{
  using namespace FourC;

  class MohrCoulombTest : public ::testing::Test
  {
   protected:
    std::shared_ptr<Mat::PlasticMohrCoulomb> make_material(double dilatancy_angle,
        double isotropic_hardening = 0.0, double saturation_cohesion = 5000.0,
        double hardening_exponent = 0.0, double youngs = 2.0e6, double cohesion = 5000.0)
    {
      Core::IO::InputParameterContainer container;
      container.add("YOUNG_MODULUS", youngs);
      container.add("POISSON_RATIO", 0.3);
      container.add("DENSITY", 2200.0);
      container.add("COHESION", cohesion);
      container.add("INTERNAL_FRICTION_ANGLE", std::numbers::pi / 6.0);
      container.add("INTERNAL_DILATANCY_ANGLE", dilatancy_angle);
      container.add("ISOTROPIC_HARDENING_MODULUS", isotropic_hardening);
      container.add("EXPONENTIAL_SATURATION_COHESION", saturation_cohesion);
      container.add("HARDENING_EXPONENT", hardening_exponent);
      container.add("TOLERANCE", 1.0e-12);
      container.add("MAX_ITERATIONS", 200);

      const int material_id = next_material_id_++;
      auto parameter = std::shared_ptr(
          Mat::make_parameter(material_id, Core::Materials::m_plmohrcoulomb, container));
      Global::Problem::instance()->materials()->set_read_from_problem(0);
      Global::Problem::instance()->materials()->insert(material_id, parameter);
      auto material = std::make_shared<Mat::PlasticMohrCoulomb>(
          dynamic_cast<Mat::PAR::PlasticMohrCoulomb*>(parameter.get()));
      parameters_.push_back(std::move(parameter));
      material->setup(1, {}, {});
      return material;
    }

    static Core::LinAlg::Tensor<double, 3, 3> diagonal_deformation_gradient(
        const std::array<double, 3>& logarithmic_strain)
    {
      Core::LinAlg::Tensor<double, 3, 3> deformation_gradient{};
      for (int i = 0; i < 3; ++i) deformation_gradient(i, i) = std::exp(logarithmic_strain[i]);
      return deformation_gradient;
    }

    static Core::LinAlg::SymmetricTensor<double, 3, 3> green_lagrange_strain(
        const Core::LinAlg::Tensor<double, 3, 3>& deformation_gradient)
    {
      return 0.5 * (Core::LinAlg::assume_symmetry(
                        Core::LinAlg::transpose(deformation_gradient) * deformation_gradient) -
                       Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
    }

    static Core::LinAlg::SymmetricTensor<double, 3, 3> evaluate(Mat::PlasticMohrCoulomb& material,
        const Core::LinAlg::Tensor<double, 3, 3>& deformation_gradient,
        Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>* tangent = nullptr)
    {
      Teuchos::ParameterList parameters;
      double total_time = 0.0;
      double time_step_size = 1.0;
      Mat::EvaluationContext<3> context{.total_time = &total_time,
          .time_step_size = &time_step_size,
          .xi = {},
          .ref_coords = nullptr};
      Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
      Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> local_tangent{};
      material.evaluate(&deformation_gradient, green_lagrange_strain(deformation_gradient),
          parameters, context, stress, local_tangent, 0, 0);
      if (tangent != nullptr) *tangent = local_tangent;
      return stress;
    }

    static Core::LinAlg::SymmetricTensor<double, 3, 3> kirchhoff_stress(
        const Core::LinAlg::Tensor<double, 3, 3>& deformation_gradient,
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& second_piola_kirchhoff_stress)
    {
      return Core::LinAlg::assume_symmetry(deformation_gradient * second_piola_kirchhoff_stress *
                                           Core::LinAlg::transpose(deformation_gradient));
    }

    int next_material_id_{1};
    std::vector<std::shared_ptr<Core::Mat::PAR::Parameter>> parameters_;
    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard_;
  };

  TEST_F(MohrCoulombTest, MatchesKratosEdgeReturn)
  {
    auto material = make_material(0.436332312998582);
    const auto deformation_gradient = diagonal_deformation_gradient({0.055, -0.122, 0.083});
    const auto tau =
        kirchhoff_stress(deformation_gradient, evaluate(*material, deformation_gradient));

    EXPECT_NEAR(tau(0, 0), -5.31271630223e4, 1.0e-5);
    EXPECT_NEAR(tau(1, 1), -1.76701997143e5, 1.0e-5);
    EXPECT_NEAR(tau(2, 2), -5.31271630223e4, 1.0e-5);
    EXPECT_NEAR(tau(0, 1), 0.0, 1.0e-12);
  }

  TEST_F(MohrCoulombTest, ReturnsToApexWithFiniteTangent)
  {
    auto material = make_material(std::numbers::pi / 6.0);
    const auto deformation_gradient = diagonal_deformation_gradient({0.1, 0.1, 0.1});
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> tangent{};
    const auto tau =
        kirchhoff_stress(deformation_gradient, evaluate(*material, deformation_gradient, &tangent));
    const double apex = 5000.0 / std::tan(std::numbers::pi / 6.0);

    for (int i = 0; i < 3; ++i) EXPECT_NEAR(tau(i, i), apex, 1.0e-6);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l) EXPECT_TRUE(std::isfinite(tangent(i, j, k, l)));
  }

  TEST_F(MohrCoulombTest, NonAssociatedFlowChangesStressReturn)
  {
    auto associated = make_material(std::numbers::pi / 6.0);
    const auto deformation_gradient = diagonal_deformation_gradient({0.055, -0.122, 0.083});
    const auto associated_tau =
        kirchhoff_stress(deformation_gradient, evaluate(*associated, deformation_gradient));

    auto non_associated = make_material(0.436332312998582);
    const auto non_associated_tau =
        kirchhoff_stress(deformation_gradient, evaluate(*non_associated, deformation_gradient));

    EXPECT_GT(std::abs(associated_tau(1, 1) - non_associated_tau(1, 1)), 1.0);
  }

  TEST_F(MohrCoulombTest, LinearAndVoceHardeningIncreaseStrength)
  {
    const auto first_deformation_gradient = diagonal_deformation_gradient({0.04, -0.08, 0.06});
    const auto second_deformation_gradient = diagonal_deformation_gradient({0.08, -0.16, 0.12});

    auto perfect = make_material(0.436332312998582);
    evaluate(*perfect, first_deformation_gradient);
    perfect->update();
    const auto perfect_tau = kirchhoff_stress(
        second_deformation_gradient, evaluate(*perfect, second_deformation_gradient));

    auto hardening = make_material(0.436332312998582, 2.0e4, 1.0e4, 20.0);
    evaluate(*hardening, first_deformation_gradient);
    hardening->update();
    const auto hardening_tau = kirchhoff_stress(
        second_deformation_gradient, evaluate(*hardening, second_deformation_gradient));

    EXPECT_GT(trace(hardening_tau), trace(perfect_tau));
  }

  TEST_F(MohrCoulombTest, HardeningCrossesExtensionEdgeApexBoundaryContinuously)
  {
    auto material = make_material(0.436332312998582, 0.0, 500.0, 10.0, 1.0e5, 100.0);

    const auto first_deformation_gradient =
        diagonal_deformation_gradient({0.011931193261700643, 0.0, 0.0});
    evaluate(*material, first_deformation_gradient);
    material->update();

    const auto second_deformation_gradient = diagonal_deformation_gradient(
        {0.01322978806459057, -0.0001121256868583982, -0.0004568607626437387});
    EXPECT_NO_THROW(evaluate(*material, second_deformation_gradient));

    Core::LinAlg::SerialDenseMatrix equivalent_plastic_strain(1, 1);
    ASSERT_TRUE(material->evaluate_output_data(
        "accumulated_equivalent_plastic_strain", equivalent_plastic_strain));
    EXPECT_NEAR(equivalent_plastic_strain(0, 0), 0.008973876448757178, 1.0e-10);
  }

  TEST_F(MohrCoulombTest, ProvidesRequestedHistoryOutputs)
  {
    auto material = make_material(0.436332312998582);
    const auto deformation_gradient = diagonal_deformation_gradient({0.055, -0.122, 0.083});
    evaluate(*material, deformation_gradient);

    Core::LinAlg::SerialDenseMatrix plastic_strain(1, 6);
    Core::LinAlg::SerialDenseMatrix equivalent_plastic_strain(1, 1);
    Core::LinAlg::SerialDenseMatrix volumetric_plastic_strain(1, 1);
    EXPECT_TRUE(material->evaluate_output_data("plastic_strain", plastic_strain));
    EXPECT_TRUE(material->evaluate_output_data(
        "accumulated_equivalent_plastic_strain", equivalent_plastic_strain));
    EXPECT_TRUE(material->evaluate_output_data(
        "accumulated_plastic_volumetric_strain", volumetric_plastic_strain));
    const double mean_plastic_strain =
        (plastic_strain(0, 0) + plastic_strain(0, 1) + plastic_strain(0, 2)) / 3.0;
    double deviatoric_plastic_strain_norm_squared = 0.0;
    for (int component = 0; component < 3; ++component)
    {
      const double deviatoric_component = plastic_strain(0, component) - mean_plastic_strain;
      deviatoric_plastic_strain_norm_squared += deviatoric_component * deviatoric_component;
    }
    const double expected_equivalent_plastic_strain =
        std::sqrt(2.0 / 3.0 * deviatoric_plastic_strain_norm_squared);
    EXPECT_NEAR(equivalent_plastic_strain(0, 0), expected_equivalent_plastic_strain, 1.0e-12);
    EXPECT_NE(volumetric_plastic_strain(0, 0), 0.0);
  }

  TEST_F(MohrCoulombTest, PacksAndUnpacksHistory)
  {
    auto material = make_material(0.436332312998582);
    const auto deformation_gradient = diagonal_deformation_gradient({0.055, -0.122, 0.083});
    const auto reference_stress = evaluate(*material, deformation_gradient);
    material->update();

    Core::Communication::PackBuffer pack_buffer;
    material->pack(pack_buffer);
    std::vector<char> packed_data;
    swap(packed_data, pack_buffer());
    Core::Communication::UnpackBuffer unpack_buffer(packed_data);
    Mat::PlasticMohrCoulomb unpacked;
    unpacked.unpack(unpack_buffer);

    FOUR_C_EXPECT_NEAR(evaluate(unpacked, deformation_gradient), reference_stress, 1.0e-8);
  }
}  // namespace
