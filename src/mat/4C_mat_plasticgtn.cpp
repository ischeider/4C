// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_plasticgtn.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_solver.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_structure_new_input.hpp"
#include "4C_utils_enum.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_local_newton.hpp"


FOUR_C_NAMESPACE_OPEN

Mat::PAR::PlasticGTN::PlasticGTN(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      youngs_(matdata.parameters.get<double>("YOUNG")),
      poissonratio_(matdata.parameters.get<double>("NUE")),
      density_(matdata.parameters.get<double>("DENS")),
      tol_(matdata.parameters.get<double>("TOL")),
      itermax_(matdata.parameters.get<int>("MAXITER")),
      yield_(matdata.parameters.get<double>("YIELD")),
      isohard_(matdata.parameters.get<double>("ISOHARD")),
      functionID_hardening_(matdata.parameters.get<int>("HARDENING_FUNC")),
      fc_(matdata.parameters.get<double>("FC")),
      kappa_(matdata.parameters.get<double>("KAPPA")),
      ef_(matdata.parameters.get<double>("EF")),
      f0_(matdata.parameters.get<double>("F0")),
      fn_(matdata.parameters.get<double>("FN")),
      sn_(matdata.parameters.get<double>("SN")),
      en_(matdata.parameters.get<double>("EN")),
      k1_(matdata.parameters.get<double>("K1")),
      k2_(matdata.parameters.get<double>("K2")),
      k3_(matdata.parameters.get<double>("K3"))
{
  if (yield_ == 0 && functionID_hardening_ == 0)
    FOUR_C_THROW(
        "You have to provide either a parameter for "
        "HARDENING_FUNC or YIELD in MAT_Struct_PlasticGTN");
}

std::shared_ptr<Core::Mat::Material> Mat::PAR::PlasticGTN::create_material()
{
  return std::make_shared<Mat::PlasticGTN>(this);
}
Mat::PlasticGTNType Mat::PlasticGTNType::instance_;

Core::Communication::ParObject* Mat::PlasticGTNType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::PlasticGTN* plastic = new Mat::PlasticGTN();
  plastic->unpack(buffer);
  return plastic;
}

Mat::PlasticGTN::PlasticGTN() : params_(nullptr) {}

Mat::PlasticGTN::PlasticGTN(Mat::PAR::PlasticGTN* params) : params_(params) {}

void Mat::PlasticGTN::pack(Core::Communication::PackBuffer& data) const
{
  int type = unique_par_object_id();
  add_to_pack(data, type);
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();
  add_to_pack(data, matid);

  // insert last converged states
  add_to_pack(data, elastic_strain_n_);
  add_to_pack(data, strain_n_);
  add_to_pack(data, stress_n_);
  add_to_pack(data, f_n_);
  add_to_pack(data, epbar_n_);

  // insert current iteration states
  add_to_pack(data, elastic_strain_n1_);
  add_to_pack(data, strain_n1_);
  add_to_pack(data, stress_n1_);
  add_to_pack(data, f_n1_);
  add_to_pack(data, epbar_n1_);
}

void Mat::PlasticGTN::unpack(Core::Communication::UnpackBuffer& buffer)
{
  isinit_ = true;

  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != nullptr)
  {
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::PlasticGTN*>(mat);
      else
        FOUR_C_THROW("Type of parameter material {} does not fit to calling type {}", mat->type(),
            material_type());

      // Extract the function for hardening again for unpack.
      const int functionID_hardening =
          params_->functionID_hardening_;  // function number for isotropic hardening
      if (functionID_hardening != 0)
      {
        hardening_function_ =
            &Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfAnything>(
                functionID_hardening);
      }
    }

    // last converged states are unpacked
    extract_from_pack(buffer, elastic_strain_n_);
    extract_from_pack(buffer, strain_n_);
    extract_from_pack(buffer, stress_n_);
    extract_from_pack(buffer, f_n_);
    extract_from_pack(buffer, epbar_n_);

    // current iteration states are unpacked
    extract_from_pack(buffer, elastic_strain_n1_);
    extract_from_pack(buffer, strain_n1_);
    extract_from_pack(buffer, stress_n1_);
    extract_from_pack(buffer, f_n1_);
    extract_from_pack(buffer, epbar_n1_);

    isinit_ = !elastic_strain_n_.empty();
  }
}

void Mat::PlasticGTN::setup(int numgp, const Discret::Elements::Fibers& fibers,
    const std::optional<Discret::Elements::CoordinateSystem>& coord_system)
{
  elastic_strain_n_.resize(numgp);
  elastic_strain_n1_.resize(numgp);
  stress_n_.resize(numgp);
  stress_n1_.resize(numgp);
  strain_n_.resize(numgp);
  strain_n1_.resize(numgp);
  f_n_.resize(numgp);
  f_n1_.resize(numgp);
  epbar_n_.resize(numgp);
  epbar_n1_.resize(numgp);

  // Extract the function for hardening only once because this is expensive.
  const int functionID_hardening =
      params_->functionID_hardening_;  // function number for isotropic hardening
  if (functionID_hardening != 0)
  {
    hardening_function_ =
        &Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfAnything>(
            functionID_hardening);
  }

  const double f0 = params_->f0_;
  std::fill(f_n_.begin(), f_n_.end(), f0);
  std::fill(f_n1_.begin(), f_n1_.end(), f0);

  isinit_ = true;
}

void Mat::PlasticGTN::update()
{
  elastic_strain_n_ = elastic_strain_n1_;
  stress_n_ = stress_n1_;
  strain_n_ = strain_n1_;
  f_n_ = f_n1_;
  epbar_n_ = epbar_n1_;
}

void Mat::PlasticGTN::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size["accumulated_plastic_strain"] = 1;
  names_and_size["damage"] = 1;
  names_and_size["elastic_strain"] = 6;
  names_and_size["plastic_strain"] = 6;
}

void Mat::PlasticGTN::evaluate(const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
    const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID)
{
  const double TOL = params_->tol_;

  // save the strain
  auto& strain_n1 = strain_n1_.at(gp);
  strain_n1 = glstrain;

  // compute incremental strain
  const Core::LinAlg::SymmetricTensor<double, 3, 3> incremental_strain =
      strain_n1 - strain_n_.at(gp);

  // compute the elastic tensor
  const double E = params_->youngs_;
  const double NU = params_->poissonratio_;
  const double G = 0.5 * E / (1 + NU);
  const double K = 1.0 / 3 * E / (1 - 2 * NU);
  const double lambda = K - 2.0 / 3.0 * G;
  constexpr auto eye = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;
  constexpr auto symmetric_identity =
      Core::LinAlg::TensorGenerators::symmetric_identity<double, 3, 3, 3, 3>;
  constexpr auto deviatoric_projection =
      symmetric_identity - 1.0 / 3.0 * Core::LinAlg::dyadic(eye, eye);
  const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> Ce =
      2.0 * G * symmetric_identity + lambda * Core::LinAlg::dyadic(eye, eye);

  // compute trial stress
  const Core::LinAlg::SymmetricTensor<double, 3, 3> stress_trial =
      stress_n_.at(gp) + Core::LinAlg::ddot(Ce, incremental_strain);
  const double p_trial = Core::LinAlg::trace(stress_trial) / 3.0;
  const Core::LinAlg::SymmetricTensor<double, 3, 3> s_trial = stress_trial - p_trial * eye;
  const double s_trial_norm = std::sqrt(Core::LinAlg::ddot(s_trial, s_trial));
  const double q_trial = std::sqrt(1.5) * s_trial_norm;

  // check the yield condition
  const double epbar_n = epbar_n_.at(gp);
  const double f_n = f_n_.at(gp);
  const double R_n = this->yield_value(epbar_n);
  const double fstar_n = this->compute_fstar(f_n);
  const double Phi = this->compute_phi(p_trial, q_trial, R_n, fstar_n);
  const double yield_cond = R_n == 0.0 ? Phi : Phi / std::abs(R_n);
  const bool is_yielded = yield_cond > TOL;

  /** plastic integration **/
  auto& stress_n1 = stress_n1_.at(gp);
  auto& f_n1 = f_n1_.at(gp);
  auto& epbar_n1 = epbar_n1_.at(gp);
  auto& elastic_strain_n1 = elastic_strain_n1_.at(gp);

  if (!is_yielded)  // elastic state
  {
    stress_n1 = stress_trial;
    f_n1 = f_n;
    epbar_n1 = epbar_n;
    elastic_strain_n1 = elastic_strain_n_.at(gp) + incremental_strain;
    stress = stress_n1;
    cmat = Ce;
  }
  else  // plastic state
  {
    // local iteration to determine p_n1, q_n1, sigmastar, f and dlambda
    constexpr int max_local_iters = 30;
    auto local_system_evaluator = [this, p_trial, q_trial, f_n, epbar_n](
                                      const Core::LinAlg::Matrix<5, 1>& x)
        -> std::pair<Core::LinAlg::Matrix<5, 1>, Core::LinAlg::Matrix<5, 5>>
    {
      const double p_n1 = x(0);
      const double q_n1 = x(1);
      const double sigmastar = x(2);
      const double f_n1 = x(3);
      const double dlambda = x(4);

      Core::LinAlg::Matrix<5, 1> h;
      this->compute_local_system(
          h, p_n1, q_n1, p_trial, q_trial, sigmastar, f_n1, f_n, dlambda, epbar_n);

      Core::LinAlg::Matrix<5, 5> dh;
      this->compute_jacobian(dh, p_n1, q_n1, sigmastar, f_n1, dlambda, epbar_n);

      return {h, dh};
    };

    Core::LinAlg::Matrix<5, 1> initial_guess;

    initial_guess(0) = p_trial;
    initial_guess(1) = q_trial;
    initial_guess(2) = R_n;
    initial_guess(3) = f_n;
    initial_guess(4) = 0.0;

    auto [x, dh] = Core::Utils::solve_local_newton_and_return_jacobian(
        local_system_evaluator, initial_guess, TOL * R_n, max_local_iters);

    const double p_n1 = x(0);
    const double q_n1 = x(1);
    const double sigmastar = x(2);
    f_n1 = x(3);
    const double dlambda = x(4);
    epbar_n1 = epbar_n + dlambda;

    // update the stress
    const Core::LinAlg::SymmetricTensor<double, 3, 3> nhat = s_trial / s_trial_norm;
    stress_n1 = sqrt(2.0 / 3) * q_n1 * nhat + p_n1 * eye;
    stress = stress_n1;

    // update elastic strain
    const double fstar = this->compute_fstar(f_n1);
    const double dphi_dsigmastar = this->compute_dphi_dsigmastar(p_n1, q_n1, sigmastar, fstar);
    const double dsigmastar_dp =
        -this->compute_dphi_dp(p_n1, q_n1, sigmastar, fstar) / dphi_dsigmastar;
    const double dsigmastar_dq =
        -this->compute_dphi_dq(p_n1, q_n1, sigmastar, fstar) / dphi_dsigmastar;
    elastic_strain_n1 = elastic_strain_n_.at(gp) + incremental_strain -
                        sqrt(1.5) * dlambda * (1.0 - f_n1) * dsigmastar_dq * nhat -
                        dlambda * (1.0 - f_n1) * dsigmastar_dp / 3.0 * eye;

    /* update the (consistent) tangent */
    Core::LinAlg::Matrix<5, 5> inv_dh;
    Core::LinAlg::Matrix<5, 1> h;

    Core::LinAlg::SymmetricTensor<double, 3, 3> dp_de{};
    Core::LinAlg::SymmetricTensor<double, 3, 3> dq_de{};
    Core::LinAlg::SymmetricTensor<double, 3, 3> dsigmastar_de{};
    Core::LinAlg::SymmetricTensor<double, 3, 3> ddlambda_de{};
    Core::LinAlg::SymmetricTensor<double, 3, 3> df_de{};

    // compute the linearization of the trial stresses w.r.t strain
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dptrial_de = K * eye;
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dqtrial_de = std::sqrt(6.0) * G * nhat;

    // compute the Jacobian and its inverse
    Core::LinAlg::FixedSizeSerialDenseSolver<5, 5> serial_dense_solver;
    inv_dh.set_copy(dh);
    serial_dense_solver.set_matrix(inv_dh);
    serial_dense_solver.invert();

    // compute the linearization of p, q, sigmastar, dlambda and f w.r.t strain
    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = i; j < 3; ++j)
      {
        h.clear();
        h(0) = dptrial_de(i, j);
        h(1) = dqtrial_de(i, j);

        x.multiply_nn(inv_dh, h);

        dp_de(i, j) = x(0);
        dq_de(i, j) = x(1);
        dsigmastar_de(i, j) = x(2);
        df_de(i, j) = x(3);
        ddlambda_de(i, j) = x(4);
      }
    }

    const double dphi_dq = this->compute_dphi_dq(p_n1, q_n1, sigmastar, fstar);

    const double d2phi_dq2 = this->compute_d2phi_dq2(p_n1, q_n1, sigmastar, fstar);
    const double d2phi_dsigmastar2 = this->compute_d2phi_dsigmastar2(p_n1, q_n1, sigmastar, fstar);
    const double d2phi_dsigmastardp =
        this->compute_d2phi_dsigmastar_dp(p_n1, q_n1, sigmastar, fstar);
    const double d2phi_dsigmastardq =
        this->compute_d2phi_dsigmastar_dq(p_n1, q_n1, sigmastar, fstar);
    const double d2phi_dsigmastardf =
        this->compute_d2phi_dsigmastar_df(p_n1, q_n1, sigmastar, f_n1);
    const double d2phi_dpdq = this->compute_d2phi_dp_dq(p_n1, q_n1, sigmastar, fstar);
    const double d2phi_dqdf = this->compute_d2phi_dq_df(p_n1, q_n1, sigmastar, f_n1);

    const double aux = dphi_dsigmastar * dphi_dsigmastar;
    const double d2sigmastar_dq2 =
        -(d2phi_dq2 * dphi_dsigmastar - d2phi_dsigmastardq * dphi_dq) / aux;
    const double d2sigmastar_dqdf =
        -(d2phi_dqdf * dphi_dsigmastar - d2phi_dsigmastardf * dphi_dq) / aux;
    const double d2sigmastar_dqdp =
        -(d2phi_dpdq * dphi_dsigmastar - d2phi_dsigmastardp * dphi_dq) / aux;
    const double d2sigmastar_dqdsigmastar =
        -(d2phi_dsigmastardq * dphi_dsigmastar - d2phi_dsigmastar2 * dphi_dq) / aux;

    // populate the consistent tangent operator
    cmat = 2 * G * (1.0 - 3 * G * dlambda / q_trial * (1.0 - f_n1) * dsigmastar_dq) *
               deviatoric_projection +
           6 * G * G * dlambda / q_trial * (1.0 - f_n1) * dsigmastar_dq *
               Core::LinAlg::dyadic(nhat, nhat) +
           sqrt(6.0) * G * dlambda * dsigmastar_dq * Core::LinAlg::dyadic(nhat, df_de) -
           sqrt(6.0) * G * (1.0 - f_n1) * dlambda * d2sigmastar_dq2 *
               Core::LinAlg::dyadic(nhat, dq_de) -
           sqrt(6.0) * G * (1.0 - f_n1) * dlambda * d2sigmastar_dqdp *
               Core::LinAlg::dyadic(nhat, dp_de) -
           sqrt(6.0) * G * (1.0 - f_n1) * dlambda * d2sigmastar_dqdsigmastar *
               Core::LinAlg::dyadic(nhat, dsigmastar_de) -
           sqrt(6.0) * G * (1.0 - f_n1) * dlambda * d2sigmastar_dqdf *
               Core::LinAlg::dyadic(nhat, df_de) -
           sqrt(6.0) * G * (1.0 - f_n1) * dsigmastar_dq * Core::LinAlg::dyadic(nhat, ddlambda_de) +
           Core::LinAlg::dyadic(eye, dp_de);
  }
}

bool Mat::PlasticGTN::evaluate_output_data(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  const auto write_strain_like_output = [&data](const std::size_t gp, const auto& strain_tensor)
  {
    data(gp, 0) = strain_tensor(0, 0);
    data(gp, 1) = strain_tensor(1, 1);
    data(gp, 2) = strain_tensor(2, 2);
    data(gp, 3) = 2.0 * strain_tensor(0, 1);
    data(gp, 4) = 2.0 * strain_tensor(1, 2);
    data(gp, 5) = 2.0 * strain_tensor(0, 2);
  };

  if (name == "accumulated_plastic_strain")
  {
    for (std::size_t gp = 0; gp < epbar_n_.size(); ++gp)
    {
      data(gp, 0) = epbar_n1_.at(gp);
    }
    return true;
  }
  else if (name == "damage")
  {
    for (std::size_t gp = 0; gp < f_n_.size(); ++gp)
    {
      data(gp, 0) = f_n1_.at(gp);
    }
    return true;
  }
  else if (name == "elastic_strain")
  {
    for (std::size_t gp = 0; gp < elastic_strain_n_.size(); ++gp)
    {
      write_strain_like_output(gp, elastic_strain_n1_.at(gp));
    }
    return true;
  }
  else if (name == "plastic_strain")
  {
    for (std::size_t gp = 0; gp < elastic_strain_n_.size(); ++gp)
    {
      const Core::LinAlg::SymmetricTensor<double, 3, 3> plastic_strain_n1 =
          strain_n1_.at(gp) - elastic_strain_n1_.at(gp);
      write_strain_like_output(gp, plastic_strain_n1);
    }
    return true;
  }
  return false;
}

double Mat::PlasticGTN::compute_fstar(const double f) const
{
  const double fC = params_->fc_;
  const double kappa = params_->kappa_;
  const double ef = params_->ef_;

  if (f <= fC - ef)
  {
    return f;
  }
  else if (f >= fC + ef)
  {
    return fC + kappa * (f - fC);
  }
  else
  {
    return (kappa - 1.0) / (4 * ef) * pow(f - fC, 2) + (kappa + 1.0) / 2 * (f - fC) + fC +
           (kappa - 1) * ef / 4;
  }
}

double Mat::PlasticGTN::compute_dfstar_df(const double f) const
{
  const double fC = params_->fc_;
  const double kappa = params_->kappa_;
  const double ef = params_->ef_;

  if (f <= fC - ef)
  {
    return 1.0;
  }
  else if (f >= fC + ef)
  {
    return kappa;
  }
  else
  {
    return (kappa - 1.0) / (2 * ef) * (f - fC) + (kappa + 1.0) / 2;
  }
}

double Mat::PlasticGTN::compute_d2fstar_df2(const double f) const
{
  const double fC = params_->fc_;
  const double kappa = params_->kappa_;
  const double ef = params_->ef_;

  if (f <= fC - ef)
  {
    return 0.0;
  }
  else if (f >= fC + ef)
  {
    return 0.0;
  }
  else
  {
    return (kappa - 1.0) / (2 * ef);
  }
}

double Mat::PlasticGTN::compute_damage_nucleation(const double alpha) const
{
  const double fN = params_->fn_;
  const double sN = params_->sn_;
  const double eN = params_->en_;
  constexpr double sqrt_2pi = 2.506628274631000502415765284811;

  return fN / (sN * sqrt_2pi) * exp(-0.5 * ((alpha - eN) / sN) * ((alpha - eN) / sN));
}

double Mat::PlasticGTN::compute_damage_nucleation_derivative(const double alpha) const
{
  const double sN = params_->sn_;
  const double eN = params_->en_;

  return -this->compute_damage_nucleation(alpha) * (alpha - eN) / (sN * sN);
}

double Mat::PlasticGTN::compute_phi(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;
  const double k3 = params_->k3_;

  return q * q / (sigmastar * sigmastar) + 2 * k1 * fstar * cosh(1.5 * k2 * p / sigmastar) - 1 -
         k3 * fstar * fstar;
}

double Mat::PlasticGTN::compute_dphi_dp(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  return 3 * k1 * k2 * fstar / sigmastar * sinh(1.5 * k2 * p / sigmastar);
}

double Mat::PlasticGTN::compute_dphi_dq(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  return 2 * q / (sigmastar * sigmastar);
}

double Mat::PlasticGTN::compute_dphi_df(
    const double p, const double q, const double sigmastar, const double f) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;
  const double k3 = params_->k3_;

  const double fstar = this->compute_fstar(f);
  const double dfstar_df = this->compute_dfstar_df(f);
  return 2 * k1 * cosh(1.5 * k2 * p / sigmastar) * dfstar_df - 2 * k3 * fstar * dfstar_df;
}

double Mat::PlasticGTN::compute_dphi_dsigmastar(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  return -2 * q * q / (sigmastar * sigmastar * sigmastar) -
         3 * k1 * k2 * fstar * p / (sigmastar * sigmastar) * sinh(1.5 * k2 * p / sigmastar);
}

double Mat::PlasticGTN::compute_d2phi_dsigmastar2(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  return 6 * pow(q, 2) / pow(sigmastar, 4) +
         6 * k1 * k2 * fstar * p / pow(sigmastar, 3) * sinh(1.5 * k2 * p / sigmastar) +
         4.5 * k1 * k2 * k2 * fstar * pow(p, 2) / pow(sigmastar, 4) *
             cosh(1.5 * k2 * p / sigmastar);
}

double Mat::PlasticGTN::compute_d2phi_dsigmastar_dp(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  return -3 * k1 * k2 * fstar / pow(sigmastar, 2) * sinh(1.5 * k2 * p / sigmastar) -
         4.5 * k1 * k2 * k2 * fstar * p / pow(sigmastar, 3) * cosh(1.5 * k2 * p / sigmastar);
}

double Mat::PlasticGTN::compute_d2phi_dsigmastar_dq(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  return -4 * q / (sigmastar * sigmastar * sigmastar);
}

double Mat::PlasticGTN::compute_d2phi_dsigmastar_df(
    const double p, const double q, const double sigmastar, const double f) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  const double dfstar_df = this->compute_dfstar_df(f);
  return -3 * k1 * k2 * p / pow(sigmastar, 2) * sinh(1.5 * k2 * p / sigmastar) * dfstar_df;
}

double Mat::PlasticGTN::compute_d2phi_dp2(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  return 4.5 * k1 * k2 * k2 * fstar / pow(sigmastar, 2) * cosh(1.5 * k2 * p / sigmastar);
}

double Mat::PlasticGTN::compute_d2phi_dp_dq(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  return 0.0;
}

double Mat::PlasticGTN::compute_d2phi_dp_df(
    const double p, const double q, const double sigmastar, const double f) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  const double dfstar_df = this->compute_dfstar_df(f);
  return 3 * k1 * k2 / sigmastar * sinh(1.5 * k2 * p / sigmastar) * dfstar_df;
}

double Mat::PlasticGTN::compute_d2phi_dq2(
    const double p, const double q, const double sigmastar, const double fstar) const
{
  return 2 / (sigmastar * sigmastar);
}

double Mat::PlasticGTN::compute_d2phi_dq_df(
    const double p, const double q, const double sigmastar, const double f) const
{
  return 0.0;
}

double Mat::PlasticGTN::compute_d2phi_df2(
    const double p, const double q, const double sigmastar, const double f) const
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;
  const double k3 = params_->k3_;

  const double fstar = this->compute_fstar(f);
  const double dfstar_df = this->compute_dfstar_df(f);
  const double d2fstar_df2 = this->compute_d2fstar_df2(f);
  return 2 * k1 * cosh(1.5 * k2 * p / sigmastar) * d2fstar_df2 - 2 * k3 * fstar * d2fstar_df2 -
         2 * k3 * dfstar_df * dfstar_df;
}

double Mat::PlasticGTN::yield_value(const double alpha) const
{
  const double yield = params_->yield_;      // initial yield stress
  const double isohard = params_->isohard_;  // linear isotropic hardening
  if (hardening_function_ == nullptr) return yield + isohard * alpha;

  const std::vector<std::pair<std::string, double>> variables{{"epsp", alpha}};
  return hardening_function_->evaluate(variables, {}, 0);
}

double Mat::PlasticGTN::yield_derivative(const double alpha) const
{
  const double isohard = params_->isohard_;  // linear isotropic hardening
  if (hardening_function_ == nullptr) return isohard;

  const std::vector<std::pair<std::string, double>> variables{{"epsp", alpha}};
  return hardening_function_->evaluate_derivative(variables, {}, 0)[0];
}

void Mat::PlasticGTN::compute_local_system(Core::LinAlg::Matrix<5, 1>& h, const double p,
    const double q, const double p_trial, const double q_trial, const double sigmastar,
    const double f, const double f_n, const double dlambda, const double epbar_n) const
{
  const double E = params_->youngs_;
  const double NU = params_->poissonratio_;

  const double G = 0.5 * E / (1 + NU);
  const double K = 1.0 / 3 * E / (1 - 2 * NU);

  const double fstar = this->compute_fstar(f);
  const double dphi_dsigmastar = this->compute_dphi_dsigmastar(p, q, sigmastar, fstar);
  const double dsigmastar_dp = -this->compute_dphi_dp(p, q, sigmastar, fstar) / dphi_dsigmastar;
  const double dsigmastar_dq = -this->compute_dphi_dq(p, q, sigmastar, fstar) / dphi_dsigmastar;
  const double A = this->compute_damage_nucleation(epbar_n + dlambda);
  h(0) = p - p_trial + (1.0 - f) * K * dsigmastar_dp * dlambda;
  h(1) = q - q_trial + 3 * (1.0 - f) * G * dsigmastar_dq * dlambda;
  h(2) = sigmastar - this->yield_value(epbar_n + dlambda);
  h(3) = this->compute_phi(p, q, sigmastar, fstar);
  h(4) = f - f_n - ((1.0 - f) * (1.0 - f) * dsigmastar_dp + A) * dlambda;
}

void Mat::PlasticGTN::compute_jacobian(Core::LinAlg::Matrix<5, 5>& jac, const double p,
    const double q, const double sigmastar, const double f, const double dlambda,
    const double epbar_n) const
{
  const double E = params_->youngs_;
  const double NU = params_->poissonratio_;

  const double G = 0.5 * E / (1 + NU);
  const double K = 1.0 / 3 * E / (1 - 2 * NU);

  const double fstar = this->compute_fstar(f);
  const double dphi_dsigmastar = this->compute_dphi_dsigmastar(p, q, sigmastar, fstar);
  const double dphi_dp = this->compute_dphi_dp(p, q, sigmastar, fstar);
  const double dphi_dq = this->compute_dphi_dq(p, q, sigmastar, fstar);
  const double dphi_df = this->compute_dphi_df(p, q, sigmastar, f);

  const double dsigmastar_dp = -dphi_dp / dphi_dsigmastar;
  const double dsigmastar_dq = -dphi_dq / dphi_dsigmastar;

  const double d2phi_dp2 = this->compute_d2phi_dp2(p, q, sigmastar, fstar);
  const double d2phi_dpdq = this->compute_d2phi_dp_dq(p, q, sigmastar, fstar);
  const double d2phi_dsigmastardp = this->compute_d2phi_dsigmastar_dp(p, q, sigmastar, fstar);
  const double d2phi_dpdf = this->compute_d2phi_dp_df(p, q, sigmastar, f);

  const double d2phi_dq2 = this->compute_d2phi_dq2(p, q, sigmastar, fstar);
  const double d2phi_dsigmastardq = this->compute_d2phi_dsigmastar_dq(p, q, sigmastar, fstar);
  const double d2phi_dqdf = this->compute_d2phi_dq_df(p, q, sigmastar, f);

  const double d2phi_dsigmastar2 = this->compute_d2phi_dsigmastar2(p, q, sigmastar, fstar);
  const double d2phi_dsigmastardf = this->compute_d2phi_dsigmastar_df(p, q, sigmastar, f);

  const double aux = dphi_dsigmastar * dphi_dsigmastar;

  const double d2sigmastar_dp2 =
      -(d2phi_dp2 * dphi_dsigmastar - d2phi_dsigmastardp * dphi_dp) / aux;
  const double d2sigmastar_dpdq =
      -(d2phi_dpdq * dphi_dsigmastar - d2phi_dsigmastardq * dphi_dp) / aux;
  const double d2sigmastar_dpdsigmastar =
      -(d2phi_dsigmastardp * dphi_dsigmastar - d2phi_dsigmastar2 * dphi_dp) / aux;
  const double d2sigmastar_dpdf =
      -(d2phi_dpdf * dphi_dsigmastar - d2phi_dsigmastardf * dphi_dp) / aux;

  const double d2sigmastar_dq2 =
      -(d2phi_dq2 * dphi_dsigmastar - d2phi_dsigmastardq * dphi_dq) / aux;
  const double d2sigmastar_dqdp =
      -(d2phi_dpdq * dphi_dsigmastar - d2phi_dsigmastardp * dphi_dq) / aux;
  const double d2sigmastar_dqdsigmastar =
      -(d2phi_dsigmastardq * dphi_dsigmastar - d2phi_dsigmastar2 * dphi_dq) / aux;
  const double d2sigmastar_dqdf =
      -(d2phi_dqdf * dphi_dsigmastar - d2phi_dsigmastardf * dphi_dq) / aux;

  jac(0, 0) = 1.0 + (1.0 - f) * K * d2sigmastar_dp2 * dlambda;
  jac(0, 1) = (1.0 - f) * K * d2sigmastar_dpdq * dlambda;
  jac(0, 2) = (1.0 - f) * K * d2sigmastar_dpdsigmastar * dlambda;
  jac(0, 3) = K * dlambda * (-dsigmastar_dp + (1.0 - f) * d2sigmastar_dpdf);
  jac(0, 4) = (1.0 - f) * K * dsigmastar_dp;

  jac(1, 0) = 3.0 * G * (1.0 - f) * d2sigmastar_dqdp * dlambda;
  jac(1, 1) = 1.0 + 3.0 * G * (1.0 - f) * d2sigmastar_dq2 * dlambda;
  jac(1, 2) = 3.0 * G * (1.0 - f) * d2sigmastar_dqdsigmastar * dlambda;
  jac(1, 3) = 3.0 * G * dlambda * (-dsigmastar_dq + (1.0 - f) * d2sigmastar_dqdf);
  jac(1, 4) = 3.0 * G * (1.0 - f) * dsigmastar_dq;

  jac(2, 0) = 0.0;
  jac(2, 1) = 0.0;
  jac(2, 2) = 1.0;
  jac(2, 3) = 0.0;
  jac(2, 4) = -this->yield_derivative(epbar_n + dlambda);

  jac(3, 0) = dphi_dp;
  jac(3, 1) = dphi_dq;
  jac(3, 2) = dphi_dsigmastar;
  jac(3, 3) = dphi_df;
  jac(3, 4) = 0.0;

  const double A = this->compute_damage_nucleation(epbar_n + dlambda);
  const double dA_dlambda = this->compute_damage_nucleation_derivative(epbar_n + dlambda);

  jac(4, 0) = -(1.0 - f) * (1.0 - f) * dlambda * d2sigmastar_dp2;
  jac(4, 1) = -(1.0 - f) * (1.0 - f) * dlambda * d2sigmastar_dpdq;
  jac(4, 2) = -(1.0 - f) * (1.0 - f) * dlambda * d2sigmastar_dpdsigmastar;
  jac(4, 3) =
      1.0 - ((1.0 - f) * (1.0 - f) * d2sigmastar_dpdf - 2.0 * (1.0 - f) * dsigmastar_dp) * dlambda;
  jac(4, 4) = -(1.0 - f) * (1.0 - f) * dsigmastar_dp - A - dA_dlambda * dlambda;
}

FOUR_C_NAMESPACE_CLOSE
