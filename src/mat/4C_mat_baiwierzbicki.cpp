// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_baiwierzbicki.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_symmetric_tensor_eigen.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_local_newton.hpp"

#include <algorithm>
#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace
{
  //! evaluate the isotropic hardening curve sigma_y(epsp)
  double evaluate_yield_stress(const Mat::PAR::BaiWierzbicki& params,
      const Core::Utils::FunctionOfAnything* hardening_function, double accplstrain)
  {
    if (hardening_function)
    {
      std::vector<std::pair<std::string, double>> dp;
      dp.emplace_back("epsp", accplstrain);
      return hardening_function->evaluate(dp, {}, 0);
    }
    return params.yield_;
  }

  //! evaluate the derivative dsigma_y/depsp of the isotropic hardening curve
  double evaluate_yield_derivative(
      const Core::Utils::FunctionOfAnything* hardening_function, double accplstrain)
  {
    if (hardening_function)
    {
      std::vector<std::pair<std::string, double>> dp;
      dp.emplace_back("epsp", accplstrain);
      return hardening_function->evaluate_derivative(dp, {}, 0)[0];
    }
    return 0.0;
  }

  //! residuum and jacobian of the return-mapping equation for the plastic multiplier.
  //! The return mapping is - by design of the MBW model - independent of the current damage,
  //! see Lian et al. (2015), eq. (37). A Perzyna overstress with viscosity visc and rate
  //! exponent rate_dependency regularizes the plastic flow, exactly as in
  //! MAT_Struct_PlasticNlnLogNeoHooke.
  std::pair<double, double> residuum_and_jacobian(const Mat::PAR::BaiWierzbicki& params,
      const Core::Utils::FunctionOfAnything* hardening_function, double Dgamma, double dt,
      double accplstrain_last, double abs_dev_KH_trial_eff)
  {
    const double ym = params.youngs_;
    const double nu = params.poissonratio_;
    const double G = ym / (2.0 * (1.0 + nu));  // undamaged shear modulus G0
    const double visc = params.visc_;
    const double eps = params.rate_dependency_;

    const double accplstrain_curr = accplstrain_last + Dgamma;
    const double y_d = evaluate_yield_stress(params, hardening_function, accplstrain_curr);
    const double dy_d_deps = evaluate_yield_derivative(hardening_function, accplstrain_curr);

    const double vfac = (visc > 0.0) ? std::pow(visc * Dgamma / dt + 1.0, eps) : 1.0;
    const double y_d_visc = y_d * vfac;
    const double dy_d_dgamma =
        (visc > 0.0) ? dy_d_deps * vfac +
                           y_d * eps * std::pow(visc * Dgamma / dt + 1.0, eps - 1.0) * visc / dt
                     : dy_d_deps;

    const double residuum =
        std::sqrt(3.0 / 2.0) * abs_dev_KH_trial_eff - 3.0 * G * Dgamma - y_d_visc;
    const double tangent = -3.0 * G - dy_d_dgamma;
    return {residuum, tangent};
  }
}  // namespace


/*----------------------------------------------------------------------*
 | constructor (public)                                                 |
 *----------------------------------------------------------------------*/
Mat::PAR::BaiWierzbicki::BaiWierzbicki(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      density_(matdata.parameters.get<double>("DENS")),
      youngs_(matdata.parameters.get<double>("YOUNG")),
      poissonratio_(matdata.parameters.get<double>("NUE")),
      yield_(matdata.parameters.get<double>("YIELD")),
      functionID_hardening_(matdata.parameters.get<int>("HARDENING_FUNC")),
      functionID_dil_(matdata.parameters.get<int>("DIL_FUNC")),
      failenergy_(matdata.parameters.get<double>("GF")),
      dcrit_(matdata.parameters.get<double>("DCRIT")),
      omega_failure_(matdata.parameters.get<double>("OMEGA_FAILURE")),
      visc_(matdata.parameters.get<double>("VISC")),
      rate_dependency_(matdata.parameters.get<double>("RATE_DEPENDENCY")),
      charlength_(matdata.parameters.get<double>("L_CHAR")),
      max_iterations_(matdata.parameters.get<int>("MAX_ITER")),
      tolerance_nr_(matdata.parameters.get<double>("TOL"))
{
  if (yield_ == 0.0 && functionID_hardening_ == 0)
    FOUR_C_THROW(
        "You have to provide either a parameter for HARDENING_FUNC or YIELD in "
        "MAT_Struct_BaiWierzbicki");

  if (functionID_dil_ != 0 && failenergy_ <= 0.0)
    FOUR_C_THROW(
        "A positive fracture energy GF is required for the damage model "
        "(i.e. whenever DIL_FUNC is given) in MAT_Struct_BaiWierzbicki");
}


/*----------------------------------------------------------------------*
 | is called in Material::Factory from read_materials()                  |
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::Mat::Material> Mat::PAR::BaiWierzbicki::create_material()
{
  return std::make_shared<Mat::BaiWierzbicki>(this);
}


Mat::BaiWierzbickiType Mat::BaiWierzbickiType::instance_;


Core::Communication::ParObject* Mat::BaiWierzbickiType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::BaiWierzbicki* bai = new Mat::BaiWierzbicki();
  bai->unpack(buffer);
  return bai;
}


/*----------------------------------------------------------------------*
 | constructors (public)                                                |
 *----------------------------------------------------------------------*/
Mat::BaiWierzbicki::BaiWierzbicki() : params_(nullptr) {}

Mat::BaiWierzbicki::BaiWierzbicki(Mat::PAR::BaiWierzbicki* params) : params_(params) {}


/*----------------------------------------------------------------------*
 | pack (public)                                                        |
 *----------------------------------------------------------------------*/
void Mat::BaiWierzbicki::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  int matid = -1;
  if (params_ != nullptr) matid = params_->id();
  add_to_pack(data, matid);

  add_to_pack(data, invplrcglast_);
  add_to_pack(data, accplstrainlast_);
  add_to_pack(data, damagelast_);
  add_to_pack(data, daminitindiclast_);
  add_to_pack(data, yieldinitlast_);
  add_to_pack(data, plastenergylast_);
  add_to_pack(data, viscenergylast_);
  add_to_pack(data, failedlast_);
  add_to_pack(data, charlength_);
}  // pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                                      |
 *----------------------------------------------------------------------*/
void Mat::BaiWierzbicki::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials()->num() != 0)
  {
    const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
    Core::Mat::PAR::Parameter* mat =
        Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
    if (mat->type() == material_type())
      params_ = static_cast<Mat::PAR::BaiWierzbicki*>(mat);
    else
      FOUR_C_THROW("Type of parameter material {} does not fit to calling type {}", mat->type(),
          material_type());

    if (params_->functionID_hardening_ != 0)
      hardening_function_ =
          &Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfAnything>(
              params_->functionID_hardening_);
    if (params_->functionID_dil_ != 0)
      dil_function_ = &Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfAnything>(
          params_->functionID_dil_);
  }

  invplrcglast_ = std::vector<Core::LinAlg::SymmetricTensor<double, 3, 3>>();
  extract_from_pack(buffer, invplrcglast_);

  accplstrainlast_ = std::vector<double>();
  extract_from_pack(buffer, accplstrainlast_);

  damagelast_ = std::vector<double>();
  extract_from_pack(buffer, damagelast_);

  daminitindiclast_ = std::vector<double>();
  extract_from_pack(buffer, daminitindiclast_);

  yieldinitlast_ = std::vector<double>();
  extract_from_pack(buffer, yieldinitlast_);

  plastenergylast_ = std::vector<double>();
  extract_from_pack(buffer, plastenergylast_);

  viscenergylast_ = std::vector<double>();
  extract_from_pack(buffer, viscenergylast_);

  failedlast_ = std::vector<double>();
  extract_from_pack(buffer, failedlast_);

  charlength_ = std::vector<double>();
  extract_from_pack(buffer, charlength_);

  const int histsize = invplrcglast_.size();
  const auto id = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;
  invplrcgcurr_.assign(histsize, id);
  accplstraincurr_.assign(histsize, 0.0);
  damagecurr_.assign(histsize, 0.0);
  daminitindiccurr_.assign(histsize, 0.0);
  yieldinitcurr_.assign(histsize, -1.0);
  plastenergycurr_.assign(histsize, 0.0);
  viscenergycurr_.assign(histsize, 0.0);
  failedcurr_.assign(histsize, 0.0);
  triaxiality_.assign(histsize, 0.0);
  lodeparam_.assign(histsize, 0.0);
  dilvalue_.assign(histsize, 0.0);
}  // unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal variables (public)                   |
 *---------------------------------------------------------------------*/
void Mat::BaiWierzbicki::setup(int numgp, const Discret::Elements::Fibers& fibers,
    const std::optional<Discret::Elements::CoordinateSystem>& coord_system)
{
  if (params_->functionID_hardening_ != 0)
    hardening_function_ =
        &Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfAnything>(
            params_->functionID_hardening_);
  if (params_->functionID_dil_ != 0)
    dil_function_ = &Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfAnything>(
        params_->functionID_dil_);

  const auto id = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;

  invplrcglast_.assign(numgp, id);
  invplrcgcurr_.assign(numgp, id);
  accplstrainlast_.assign(numgp, 0.0);
  accplstraincurr_.assign(numgp, 0.0);
  damagelast_.assign(numgp, 0.0);
  damagecurr_.assign(numgp, 0.0);
  daminitindiclast_.assign(numgp, 0.0);
  daminitindiccurr_.assign(numgp, 0.0);
  yieldinitlast_.assign(numgp, -1.0);
  yieldinitcurr_.assign(numgp, -1.0);
  plastenergylast_.assign(numgp, 0.0);
  plastenergycurr_.assign(numgp, 0.0);
  viscenergylast_.assign(numgp, 0.0);
  viscenergycurr_.assign(numgp, 0.0);
  failedlast_.assign(numgp, 0.0);
  failedcurr_.assign(numgp, 0.0);
  charlength_.assign(numgp, -1.0);
  triaxiality_.assign(numgp, 0.0);
  lodeparam_.assign(numgp, 0.0);
  dilvalue_.assign(numgp, 0.0);

  isinit_ = true;
}  // setup()


/*----------------------------------------------------------------------*
 | update internal variables                                            |
 *----------------------------------------------------------------------*/
void Mat::BaiWierzbicki::update()
{
  invplrcglast_ = invplrcgcurr_;
  accplstrainlast_ = accplstraincurr_;
  damagelast_ = damagecurr_;
  daminitindiclast_ = daminitindiccurr_;
  yieldinitlast_ = yieldinitcurr_;
  plastenergylast_ = plastenergycurr_;
  viscenergylast_ = viscenergycurr_;
  failedlast_ = failedcurr_;

  const int histsize = invplrcglast_.size();
  const auto id = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;
  invplrcgcurr_.assign(histsize, id);
  accplstraincurr_.assign(histsize, 0.0);
  damagecurr_.assign(histsize, 0.0);
  daminitindiccurr_.assign(histsize, 0.0);
  yieldinitcurr_.assign(histsize, -1.0);
  plastenergycurr_.assign(histsize, 0.0);
  viscenergycurr_.assign(histsize, 0.0);
  failedcurr_.assign(histsize, 0.0);
}  // update()


/*----------------------------------------------------------------------*
 | calculate stress and constitutive tensor                             |
 *----------------------------------------------------------------------*/
void Mat::BaiWierzbicki::evaluate(const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
    const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID)
{
  // elastic (undamaged) material data
  const double ym = params_->youngs_;
  const double nu = params_->poissonratio_;
  const double G0 = ym / (2.0 * (1.0 + nu));        // undamaged shear modulus
  const double K0 = ym / (3.0 * (1.0 - 2.0 * nu));  // undamaged bulk modulus

  // damage / regularization data
  const double Gf = params_->failenergy_;
  const double Dcr = params_->dcrit_;
  const double omega_fail = params_->omega_failure_;
  const double visc = params_->visc_;
  const double rate_dep = params_->rate_dependency_;

  FOUR_C_ASSERT(context.time_step_size, "Time step size not given in evaluation context.");
  const double dt = *context.time_step_size;

  // characteristic length: user parameter, or characteristic Gauss point size from the element
  double L = params_->charlength_;
  if (L <= 0.0) L = (context.gp_char_length != nullptr) ? *context.gp_char_length : 0.0;
  charlength_.at(gp) = L;

  const double detF = Core::LinAlg::det(*defgrad);
  const auto invdefgrd = Core::LinAlg::inv(*defgrad);

  // The logarithmic-strain framework is only valid for positive stretches / detF. During the
  // global nonlinear iteration a trial displacement iterate may momentarily invert an element,
  // producing detF <= 0 or non-positive eigenvalues of the elastic left Cauchy-Green tensor.
  // With floating-point trapping enabled, feeding such values to std::log would raise a fatal
  // SIGFPE and bypass the time-step adaptation. We therefore floor the arguments of std::log by
  // a tiny positive number: converged, physical states are unaffected (detF, eigenvalues > 0),
  // while non-physical trial states yield a bounded response that the solver rejects.
  constexpr double log_floor = 1.0e-12;
  const double log_detF = std::log(std::max(detF, log_floor));

  const auto id2 = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;
  const auto Idev = id2 - 1.0 / 3.0 * Core::LinAlg::TensorGenerators::ones<double, 3, 3>;

  // effective (undamaged) principal material tangent, initialized elastic
  Core::LinAlg::SymmetricTensor<double, 3, 3> D_ep_principal_eff =
      Core::LinAlg::TensorGenerators::ones<double, 3, 3> * K0;
  D_ep_principal_eff += 2.0 * G0 * Idev;

  // ---------------------------------------------------- trial elastic state
  // elastic left Cauchy-Green tensor Be_trial = F C_{p,n}^{-1} F^T
  const auto Be_trial = Core::LinAlg::assume_symmetry(
      *defgrad * invplrcglast_.at(gp) * Core::LinAlg::transpose(*defgrad));

  const auto& [lambda_trial_square, n] = Core::LinAlg::eig(Be_trial);

  std::vector<Core::LinAlg::Tensor<double, 3>> spatial_principal_directions(3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) (spatial_principal_directions.at(i))(j) = n(j, i);

  std::vector<Core::LinAlg::Tensor<double, 3>> material_principal_directions(3);
  for (int i = 0; i < 3; i++)
    material_principal_directions.at(i) = invdefgrd * spatial_principal_directions.at(i);

  // effective (undamaged) deviatoric Kirchhoff stress at trial state
  Core::LinAlg::Tensor<double, 3> dev_KH_trial_eff;
  for (int i = 0; i < 3; i++)
    dev_KH_trial_eff(i) =
        G0 * std::log(std::max(lambda_trial_square[i], log_floor)) - 2.0 / 3.0 * G0 * log_detF;

  // effective (undamaged) pressure (equal at trial and final state)
  const double pressure_eff = K0 * log_detF;

  const double abs_dev_eff = Core::LinAlg::norm2(dev_KH_trial_eff);

  // effective deviatoric Kirchhoff stress at final state (updated in the plastic branch)
  Core::LinAlg::Tensor<double, 3> dev_KH_eff = dev_KH_trial_eff;
  Core::LinAlg::SymmetricTensor<double, 3, 3> Be = Be_trial;

  const double D_n = damagelast_.at(gp);
  const bool already_failed = failedlast_.at(gp) > 0.5;

  double Dgamma = 0.0;
  double eta = 0.0;
  double thetabar = 0.0;
  double D_np1 = D_n;
  double w_np1 = daminitindiclast_.at(gp);
  double sigy0 = yieldinitlast_.at(gp);
  double dWp = 0.0;
  double dWv = 0.0;
  bool failed_np1 = already_failed;

  if (!already_failed)
  {
    const double y_d_n =
        evaluate_yield_stress(*params_, hardening_function_, accplstrainlast_.at(gp));
    const double f_trial = std::sqrt(3.0 / 2.0) * abs_dev_eff - y_d_n;

    if (f_trial > 0.0)  // -------------------------------------- plastic step
    {
      auto residuumAndJacobian = [&](double dg)
      {
        return residuum_and_jacobian(
            *params_, hardening_function_, dg, dt, accplstrainlast_.at(gp), abs_dev_eff);
      };
      Dgamma = Core::Utils::solve_local_newton(
          residuumAndJacobian, 0.0, params_->tolerance_nr_, params_->max_iterations_);

      // flow vector in principal space (effective)
      Core::LinAlg::Tensor<double, 3> flow_vector(dev_KH_trial_eff);
      flow_vector /= std::sqrt(2.0 / 3.0) * abs_dev_eff;

      // deviatoric stress return mapping (effective)
      dev_KH_eff =
          (1.0 - 2.0 * G0 * Dgamma / (std::sqrt(2.0 / 3.0) * abs_dev_eff)) * dev_KH_trial_eff;

      // strain return mapping
      Be = Core::LinAlg::SymmetricTensor<double, 3, 3>{};
      for (int i = 0; i < 3; i++)
        Be += lambda_trial_square[i] / std::exp(2.0 * flow_vector(i) * Dgamma) *
              Core::LinAlg::self_dyadic(spatial_principal_directions.at(i));

      const double accplstrain_np1 = accplstrainlast_.at(gp) + Dgamma;
      const double y_d_curr = evaluate_yield_stress(*params_, hardening_function_, accplstrain_np1);

      // stress-state descriptors from the converged effective deviatoric stress
      const double sigma_vm = std::sqrt(3.0 / 2.0) * Core::LinAlg::norm2(dev_KH_eff);
      if (sigma_vm > 1.0e-12)
      {
        eta = pressure_eff / sigma_vm;
        const double J3 = dev_KH_eff(0) * dev_KH_eff(1) * dev_KH_eff(2);
        double xi = 13.5 * J3 / (sigma_vm * sigma_vm * sigma_vm);
        xi = std::clamp(xi, -1.0, 1.0);
        thetabar = 1.0 - 2.0 / M_PI * std::acos(xi);
      }
      triaxiality_.at(gp) = eta;
      lodeparam_.at(gp) = thetabar;

      // damage initiation and evolution
      double dil_value = 0.0;
      if (dil_function_)
      {
        std::vector<std::pair<std::string, double>> dpar;
        dpar.emplace_back("eta", eta);
        dpar.emplace_back("theta", thetabar);
        dil_value = dil_function_->evaluate(dpar, {}, 0);
        if (dil_value < 1.0e-12) dil_value = 1.0e-12;

        const double w_last = daminitindiclast_.at(gp);
        w_np1 = w_last + Dgamma / dil_value;

        if (sigy0 < 0.0)  // damage not yet initiated
        {
          if (w_np1 >= 1.0)
          {
            sigy0 = y_d_curr;  // freeze yield stress at damage initiation
            const double frac = (w_np1 - 1.0) / std::max(w_np1 - w_last, 1.0e-12);
            D_np1 = D_n + sigy0 * L / (2.0 * Gf) * frac * Dgamma;
          }
        }
        else  // damage already initiated
        {
          D_np1 = D_n + sigy0 * L / (2.0 * Gf) * Dgamma;
        }
        if (D_np1 > Dcr) D_np1 = Dcr;
      }
      dilvalue_.at(gp) = dil_value;

      // energies (damaged, converged flow stress omega * y_d * vfac)
      const double omega_np1 = 1.0 - D_np1;
      const double vfac = (visc > 0.0) ? std::pow(visc * Dgamma / dt + 1.0, rate_dep) : 1.0;
      dWp = omega_np1 * y_d_curr * Dgamma;
      dWv = omega_np1 * y_d_curr * (vfac - 1.0) * Dgamma;

      // consistent effective plastic tangent (as in MAT_Struct_PlasticNlnLogNeoHooke)
      const auto grad = residuumAndJacobian(Dgamma);
      const double dy_d_dgamma = -grad.second - 3.0 * G0;  // = dsigma_y,visc/dDgamma
      const double fac1 = -4.0 * G0 * G0 * Dgamma / (std::sqrt(2.0 / 3.0) * abs_dev_eff);
      D_ep_principal_eff += fac1 * Idev;
      const double fac2 =
          4.0 * G0 * G0 *
          (std::sqrt(2.0 / 3.0) * Dgamma / abs_dev_eff - 1.0 / (3.0 * G0 + dy_d_dgamma));
      D_ep_principal_eff += fac2 * Core::LinAlg::self_dyadic(flow_vector);
    }
    else  // -------------------------------------------------- elastic step
    {
      const double sigma_vm = std::sqrt(3.0 / 2.0) * abs_dev_eff;
      if (sigma_vm > 1.0e-12)
      {
        eta = pressure_eff / sigma_vm;
        const double J3 = dev_KH_eff(0) * dev_KH_eff(1) * dev_KH_eff(2);
        double xi = 13.5 * J3 / (sigma_vm * sigma_vm * sigma_vm);
        xi = std::clamp(xi, -1.0, 1.0);
        thetabar = 1.0 - 2.0 / M_PI * std::acos(xi);
      }
      triaxiality_.at(gp) = eta;
      lodeparam_.at(gp) = thetabar;
      if (dil_function_)
      {
        std::vector<std::pair<std::string, double>> dpar;
        dpar.emplace_back("eta", eta);
        dpar.emplace_back("theta", thetabar);
        dilvalue_.at(gp) = dil_function_->evaluate(dpar, {}, 0);
      }
    }

    // mark Gauss point as failed once the integrity drops to the failure limit
    if ((1.0 - D_np1) <= omega_fail) failed_np1 = true;
  }
  else  // ----------------------------------------------- already failed GP
  {
    // stress is omega_fail times the elastic response; internal variables are frozen
    triaxiality_.at(gp) = triaxiality_.at(gp);
    lodeparam_.at(gp) = lodeparam_.at(gp);
  }

  // integrity used for stress/tangent, floored at the failure limit to avoid singularity
  const double omega_used = std::max(1.0 - D_np1, omega_fail);

  // ----------------------------------------------------- actual (damaged) stress
  Core::LinAlg::Tensor<double, 3> KH_actual;
  for (int i = 0; i < 3; i++) KH_actual(i) = omega_used * (dev_KH_eff(i) + pressure_eff);

  Core::LinAlg::SymmetricTensor<double, 3, 3> D_ep_principal = omega_used * D_ep_principal_eff;

  stress = Core::LinAlg::SymmetricTensor<double, 3, 3>{};
  for (int i = 0; i < 3; i++)
    stress += KH_actual(i) * Core::LinAlg::self_dyadic(material_principal_directions.at(i));

  // ------------------------------------------------------------ tangent modulus
  cmat = Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>{};
  std::vector<Core::LinAlg::Tensor<double, 3, 3>> mat_princ_direction_dyad(3);
  for (int a = 0; a < 3; a++)
    mat_princ_direction_dyad.at(a) = Core::LinAlg::dyadic(
        material_principal_directions.at(a), material_principal_directions.at(a));

  for (int a = 0; a < 3; a++)
  {
    cmat += -2.0 * KH_actual(a) *
            Core::LinAlg::assume_symmetry(Core::LinAlg::dyadic(
                mat_princ_direction_dyad.at(a), mat_princ_direction_dyad.at(a)));

    for (int b = 0; b < 3; b++)
    {
      cmat += D_ep_principal(a, b) *
              Core::LinAlg::assume_symmetry(Core::LinAlg::dyadic(
                  mat_princ_direction_dyad.at(a), mat_princ_direction_dyad.at(b)));

      if (a != b)
      {
        const double fac =
            (lambda_trial_square[a] != lambda_trial_square[b])
                ? (KH_actual(a) * lambda_trial_square[b] - KH_actual(b) * lambda_trial_square[a]) /
                      (lambda_trial_square[a] - lambda_trial_square[b])
                : 0.5 * (D_ep_principal(b, b) - D_ep_principal(a, b)) - KH_actual(b);
        const auto tmp1 = Core::LinAlg::dyadic(
            material_principal_directions.at(a), material_principal_directions.at(b));
        cmat += fac * Core::LinAlg::assume_symmetry(Core::LinAlg::dyadic(tmp1, tmp1));
        cmat += fac * Core::LinAlg::assume_symmetry(Core::LinAlg::dyadic(tmp1, transpose(tmp1)));
      }
    }
  }

  // ---------------------------------------------------------- update history
  if (already_failed)
  {
    invplrcgcurr_.at(gp) = invplrcglast_.at(gp);
    accplstraincurr_.at(gp) = accplstrainlast_.at(gp);
    damagecurr_.at(gp) = damagelast_.at(gp);
    daminitindiccurr_.at(gp) = daminitindiclast_.at(gp);
    yieldinitcurr_.at(gp) = yieldinitlast_.at(gp);
    plastenergycurr_.at(gp) = plastenergylast_.at(gp);
    viscenergycurr_.at(gp) = viscenergylast_.at(gp);
    failedcurr_.at(gp) = 1.0;
  }
  else
  {
    invplrcgcurr_.at(gp) =
        Core::LinAlg::assume_symmetry(invdefgrd * Be * Core::LinAlg::transpose(invdefgrd));
    accplstraincurr_.at(gp) = accplstrainlast_.at(gp) + Dgamma;
    damagecurr_.at(gp) = D_np1;
    daminitindiccurr_.at(gp) = w_np1;
    yieldinitcurr_.at(gp) = sigy0;
    plastenergycurr_.at(gp) = plastenergylast_.at(gp) + dWp;
    viscenergycurr_.at(gp) = viscenergylast_.at(gp) + dWv;
    failedcurr_.at(gp) = failed_np1 ? 1.0 : 0.0;
  }
}  // evaluate()


/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void Mat::BaiWierzbicki::vis_names(std::map<std::string, int>& names) const
{
  names["accumulatedstrain"] = 1;
  names["damage"] = 1;
  names["triaxiality"] = 1;
  names["lode_angle_parameter"] = 1;
  names["damage_initiation_locus"] = 1;
  names["gp_size"] = 1;
  names["plastic_dissipated_energy"] = 1;
  names["viscous_artificial_energy"] = 1;
}  // vis_names()


/*---------------------------------------------------------------------*
 | return names of visualization data for direct VTK output            |
 *---------------------------------------------------------------------*/
void Mat::BaiWierzbicki::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size["accumulated_plastic_strain"] = 1;
  names_and_size["damage"] = 1;
  names_and_size["triaxiality"] = 1;
  names_and_size["lode_angle_parameter"] = 1;
  names_and_size["damage_initiation_locus"] = 1;
  names_and_size["gp_size"] = 1;
  names_and_size["plastic_dissipated_energy"] = 1;
  names_and_size["viscous_artificial_energy"] = 1;
  names_and_size["failed_flag"] = 1;
}


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool Mat::BaiWierzbicki::vis_data(
    const std::string& name, std::vector<double>& data, int numgp, int eleID) const
{
  auto average = [&](const std::vector<double>& v)
  {
    double temp = 0.0;
    for (int gp = 0; gp < numgp; ++gp) temp += v.at(gp);
    data[0] = temp / numgp;
  };

  if (name == "accumulatedstrain")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
    average(accplstrainlast_);
    return true;
  }
  if (name == "damage")
  {
    average(damagelast_);
    return true;
  }
  if (name == "triaxiality")
  {
    average(triaxiality_);
    return true;
  }
  if (name == "lode_angle_parameter")
  {
    average(lodeparam_);
    return true;
  }
  if (name == "damage_initiation_locus")
  {
    average(dilvalue_);
    return true;
  }
  if (name == "gp_size")
  {
    average(charlength_);
    return true;
  }
  if (name == "plastic_dissipated_energy")
  {
    average(plastenergylast_);
    return true;
  }
  if (name == "viscous_artificial_energy")
  {
    average(viscenergylast_);
    return true;
  }
  return false;
}  // vis_data()


bool Mat::BaiWierzbicki::evaluate_output_data(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  auto fill = [&](const std::vector<double>& v)
  {
    for (std::size_t gp = 0; gp < v.size(); ++gp) data(gp, 0) = v.at(gp);
  };

  if (name == "accumulated_plastic_strain")
  {
    fill(accplstraincurr_);
    return true;
  }
  if (name == "damage")
  {
    fill(damagecurr_);
    return true;
  }
  if (name == "triaxiality")
  {
    fill(triaxiality_);
    return true;
  }
  if (name == "lode_angle_parameter")
  {
    fill(lodeparam_);
    return true;
  }
  if (name == "damage_initiation_locus")
  {
    fill(dilvalue_);
    return true;
  }
  if (name == "gp_size")
  {
    fill(charlength_);
    return true;
  }
  if (name == "plastic_dissipated_energy")
  {
    fill(plastenergycurr_);
    return true;
  }
  if (name == "viscous_artificial_energy")
  {
    fill(viscenergycurr_);
    return true;
  }
  if (name == "failed_flag")
  {
    fill(failedcurr_);
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
