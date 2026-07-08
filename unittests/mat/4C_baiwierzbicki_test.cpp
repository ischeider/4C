// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_comm_pack_buffer.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_tensor_conversion.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_baiwierzbicki.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <cmath>

namespace
{
  using namespace FourC;

  //! build a parameter container for a constant-yield, no-damage BaiWierzbicki material
  Core::IO::InputParameterContainer make_container(double young, double nue, double yield)
  {
    Core::IO::InputParameterContainer container;
    container.add("DENS", 1.0);
    container.add("YOUNG", young);
    container.add("NUE", nue);
    container.add("YIELD", yield);
    container.add("HARDENING_FUNC", 0);
    container.add("DIL_FUNC", 0);
    container.add("GF", 0.0);
    container.add("DCRIT", 1.0);
    container.add("OMEGA_FAILURE", 1.0e-4);
    container.add("VISC", 0.0);
    container.add("RATE_DEPENDENCY", 0.0);
    container.add("L_CHAR", 1.0);
    container.add("MAX_ITER", 50);
    container.add("TOL", 1.0e-12);
    return container;
  }

  class BaiWierzbickiTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      const Core::IO::InputParameterContainer container = make_container(1000.0, 0.3, 1.0e30);

      param_ = std::shared_ptr(
          Mat::make_parameter(1, Core::Materials::MaterialType::m_baiwierzbicki, container));

      Global::Problem& problem = (*Global::Problem::instance());
      problem.materials()->set_read_from_problem(0);
      problem.materials()->insert(1, param_);
      bai_ = std::make_shared<Mat::BaiWierzbicki>(
          dynamic_cast<Mat::PAR::BaiWierzbicki*>(param_.get()));
    }

    std::shared_ptr<Core::Mat::PAR::Parameter> param_;
    Core::Communication::PackBuffer data;
    std::shared_ptr<Mat::BaiWierzbicki> bai_;

    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard;
  };

  //! For a purely elastic small deformation the log-strain material must reproduce
  //! the linear-elastic 2nd Piola-Kirchhoff stress to high accuracy.
  TEST_F(BaiWierzbickiTest, ElasticStressMatchesLinearElasticity)
  {
    bai_->setup(1, {}, {});

    const double young = 1000.0;
    const double nue = 0.3;
    const double lambda = young * nue / ((1.0 + nue) * (1.0 - 2.0 * nue));
    const double mu = young / (2.0 * (1.0 + nue));

    // small uniaxial-in-x Green-Lagrange strain
    const double eps = 1.0e-5;
    Core::LinAlg::SymmetricTensor<double, 3, 3> glstrain{};
    glstrain(0, 0) = eps;

    const double trace = eps;
    Core::LinAlg::SymmetricTensor<double, 3, 3> ref_stress{};
    for (int i = 0; i < 3; ++i) ref_stress(i, i) = lambda * trace;
    ref_stress(0, 0) += 2.0 * mu * eps;

    Core::LinAlg::Tensor<double, 3, 3> defgrad{};
    defgrad(0, 0) = std::sqrt(1.0 + 2.0 * eps);
    defgrad(1, 1) = 1.0;
    defgrad(2, 2) = 1.0;

    Teuchos::ParameterList paras;
    double total_time = 0.0;
    double time_step_size = 1.0;
    double gp_char = 1.0;
    Mat::EvaluationContext<3> context{.total_time = &total_time,
        .time_step_size = &time_step_size,
        .xi = nullptr,
        .ref_coords = nullptr,
        .gp_char_length = &gp_char};

    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat{};
    Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
    bai_->evaluate(&defgrad, glstrain, paras, context, stress, cmat, 0, 0);

    // tolerance is loose because the log-strain framework only matches to first order in eps
    FOUR_C_EXPECT_NEAR(stress, ref_stress, 1.0e-4);
  }

  //! test pack/unpack round-trip of the material history
  TEST_F(BaiWierzbickiTest, PackUnpack)
  {
    bai_->setup(1, {}, {});

    bai_->pack(data);
    std::vector<char> dataSend;
    swap(dataSend, data());

    Mat::BaiWierzbicki restored;
    Core::Communication::UnpackBuffer buffer(dataSend);
    restored.unpack(buffer);

    EXPECT_EQ(restored.unique_par_object_id(), bai_->unique_par_object_id());
    EXPECT_EQ(restored.material_type(), Core::Materials::m_baiwierzbicki);
  }

  //! J2 return mapping: for a deviatoric (isochoric) load beyond yield, the resulting von Mises
  //! Kirchhoff stress must equal the (constant) yield stress and plasticity must be active.
  TEST(BaiWierzbickiReturnMapTest, PlasticYieldStressIsRecovered)
  {
    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard;

    const double young = 1000.0;
    const double nue = 0.3;
    const double yield = 1.0;

    Core::IO::InputParameterContainer container;
    container.add("DENS", 1.0);
    container.add("YOUNG", young);
    container.add("NUE", nue);
    container.add("YIELD", yield);
    container.add("HARDENING_FUNC", 0);
    container.add("DIL_FUNC", 0);
    container.add("GF", 0.0);
    container.add("DCRIT", 1.0);
    container.add("OMEGA_FAILURE", 1.0e-4);
    container.add("VISC", 0.0);
    container.add("RATE_DEPENDENCY", 0.0);
    container.add("L_CHAR", 1.0);
    container.add("MAX_ITER", 50);
    container.add("TOL", 1.0e-12);

    auto param =
        std::shared_ptr(Mat::make_parameter(1, Core::Materials::m_baiwierzbicki, container));
    Global::Problem& problem = (*Global::Problem::instance());
    problem.materials()->set_read_from_problem(0);
    problem.materials()->insert(1, param);
    auto bai =
        std::make_shared<Mat::BaiWierzbicki>(dynamic_cast<Mat::PAR::BaiWierzbicki*>(param.get()));
    bai->setup(1, {}, {});

    // isochoric shear deformation gradient large enough to yield
    const double gamma = 0.05;  // engineering shear
    Core::LinAlg::Tensor<double, 3, 3> defgrad{};
    defgrad(0, 0) = 1.0;
    defgrad(1, 1) = 1.0;
    defgrad(2, 2) = 1.0;
    defgrad(0, 1) = gamma;

    // Green-Lagrange strain E = 0.5 (F^T F - I)
    Core::LinAlg::SymmetricTensor<double, 3, 3> glstrain{};
    glstrain(0, 1) = 0.5 * gamma;
    glstrain(1, 1) = 0.5 * gamma * gamma;

    Teuchos::ParameterList paras;
    double total_time = 0.0;
    double time_step_size = 1.0;
    double gp_char = 1.0;
    Mat::EvaluationContext<3> context{.total_time = &total_time,
        .time_step_size = &time_step_size,
        .xi = nullptr,
        .ref_coords = nullptr,
        .gp_char_length = &gp_char};

    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat{};
    Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
    bai->evaluate(&defgrad, glstrain, paras, context, stress, cmat, 0, 0);

    // Kirchhoff stress tau = F S F^T ; von Mises of tau must equal the yield stress
    Core::LinAlg::Matrix<3, 3> F = Core::LinAlg::make_matrix(defgrad);
    Core::LinAlg::Matrix<3, 3> S = Core::LinAlg::make_matrix(Core::LinAlg::get_full(stress));
    Core::LinAlg::Matrix<3, 3> FS;
    FS.multiply_nn(F, S);
    Core::LinAlg::Matrix<3, 3> tau;
    tau.multiply_nt(FS, F);

    const double trace = tau(0, 0) + tau(1, 1) + tau(2, 2);
    Core::LinAlg::Matrix<3, 3> dev(tau);
    for (int i = 0; i < 3; ++i) dev(i, i) -= trace / 3.0;
    double j2 = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) j2 += 0.5 * dev(i, j) * dev(i, j);
    const double vm = std::sqrt(3.0 * j2);

    EXPECT_NEAR(vm, yield, 1.0e-6);
    EXPECT_GT(bai->accumulated_strain(0), 0.0);
    EXPECT_DOUBLE_EQ(bai->damage(0), 0.0);
  }
}  // namespace
