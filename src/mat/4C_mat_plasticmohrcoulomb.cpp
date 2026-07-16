// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_plasticmohrcoulomb.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_symmetric_tensor_eigen.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_base.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

namespace
{
  using Vector3 = std::array<double, 3>;
  using Matrix3 = std::array<Vector3, 3>;

  constexpr double numerical_tangent_perturbation = 1.0e-7;

  double dot(const Vector3& a, const Vector3& b)
  {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
  }

  double norm(const Vector3& a) { return std::sqrt(dot(a, a)); }

  Vector3 deviator(const Vector3& vector)
  {
    const double mean = (vector[0] + vector[1] + vector[2]) / 3.0;
    return {vector[0] - mean, vector[1] - mean, vector[2] - mean};
  }

  double equivalent_plastic_strain_increment(const Vector3& plastic_strain_increment)
  {
    const Vector3 deviatoric_increment = deviator(plastic_strain_increment);
    return std::sqrt(2.0 / 3.0 * dot(deviatoric_increment, deviatoric_increment));
  }

  Vector3 operator-(const Vector3& a, const Vector3& b)
  {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  }

  Vector3 cross(const Vector3& a, const Vector3& b)
  {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
  }

  double safe_denominator(double value)
  {
    constexpr double minimum = 1.0e-14;
    if (std::abs(value) >= minimum) return value;
    return std::copysign(minimum, value == 0.0 ? 1.0 : value);
  }

  struct PrincipalReturn
  {
    Vector3 stress{};
    Vector3 plastic_strain_increment{};
    int region{};
  };

  Matrix3 elastic_normal_matrix(double youngs, double poisson_ratio)
  {
    const double diagonal =
        youngs * (1.0 - poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    const double off_diagonal =
        youngs * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    Matrix3 matrix{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) matrix[i][j] = (i == j) ? diagonal : off_diagonal;
    return matrix;
  }

  Vector3 elastic_normal_compliance_times(
      const Vector3& stress, double youngs, double poisson_ratio)
  {
    const double trace_stress = stress[0] + stress[1] + stress[2];
    Vector3 strain{};
    for (int i = 0; i < 3; ++i)
      strain[i] = ((1.0 + poisson_ratio) * stress[i] - poisson_ratio * trace_stress) / youngs;
    return strain;
  }

  PrincipalReturn return_to_mohr_coulomb(const Vector3& trial_stress, double cohesion,
      double friction_angle, double dilatancy_angle, double youngs, double poisson_ratio)
  {
    const double friction_coefficient =
        (1.0 + std::sin(friction_angle)) / (1.0 - std::sin(friction_angle));
    const double dilatancy_coefficient =
        (1.0 + std::sin(dilatancy_angle)) / (1.0 - std::sin(dilatancy_angle));
    const double cohesion_coefficient = 2.0 * cohesion * std::sqrt(friction_coefficient);
    const double state_function =
        friction_coefficient * trial_stress[0] - trial_stress[2] - cohesion_coefficient;

    PrincipalReturn result{.stress = trial_stress};
    if (state_function <= 0.0) return result;

    const double apex = cohesion_coefficient / (friction_coefficient - 1.0);
    const Matrix3 elastic = elastic_normal_matrix(youngs, poisson_ratio);

    const double den_p = safe_denominator(
        friction_coefficient * (elastic[0][0] * dilatancy_coefficient - elastic[0][2]) -
        elastic[2][0] * dilatancy_coefficient + elastic[2][2]);
    const Vector3 return_plane = {(elastic[0][0] * dilatancy_coefficient - elastic[0][2]) / den_p,
        (elastic[1][0] * dilatancy_coefficient - elastic[1][2]) / den_p,
        (elastic[2][0] * dilatancy_coefficient - elastic[2][2]) / den_p};

    const Vector3 trial_from_apex = {
        trial_stress[0] - apex, trial_stress[1] - apex, trial_stress[2] - apex};
    const Vector3 boundary_i_ii = {return_plane[1] * friction_coefficient - return_plane[2],
        return_plane[2] - return_plane[0] * friction_coefficient,
        return_plane[0] - return_plane[1]};
    const Vector3 boundary_i_iii = {friction_coefficient * (return_plane[1] - return_plane[2]),
        return_plane[2] - return_plane[0] * friction_coefficient,
        return_plane[0] * friction_coefficient - return_plane[1]};
    const double p_i_ii = dot(boundary_i_ii, trial_from_apex);
    const double p_i_iii = dot(boundary_i_iii, trial_from_apex);

    const double den_p2 = safe_denominator(
        friction_coefficient * (elastic[1][1] * dilatancy_coefficient - elastic[1][2]) -
        elastic[1][2] * dilatancy_coefficient + elastic[2][2]);
    const Vector3 return_secondary_ii = {
        (elastic[0][1] * dilatancy_coefficient - elastic[0][2]) / den_p2,
        (elastic[1][1] * dilatancy_coefficient - elastic[1][2]) / den_p2,
        (elastic[2][1] * dilatancy_coefficient - elastic[2][2]) / den_p2};
    const Vector3 normal_ii = cross(return_plane, return_secondary_ii);
    const double t1 =
        dot(normal_ii, trial_from_apex) /
        safe_denominator(normal_ii[0] + normal_ii[1] + friction_coefficient * normal_ii[2]);

    const double den_p3 = safe_denominator(
        friction_coefficient * (elastic[0][0] * dilatancy_coefficient - elastic[0][1]) -
        elastic[1][0] * dilatancy_coefficient + elastic[1][1]);
    const Vector3 return_secondary_iii = {
        (elastic[0][0] * dilatancy_coefficient - elastic[0][1]) / den_p3,
        (elastic[1][0] * dilatancy_coefficient - elastic[1][1]) / den_p3,
        (elastic[2][0] * dilatancy_coefficient - elastic[2][1]) / den_p3};
    const Vector3 normal_iii = cross(return_plane, return_secondary_iii);
    const double t2 = dot(normal_iii, trial_from_apex) /
                      safe_denominator(normal_iii[0] + friction_coefficient * normal_iii[1] +
                                       friction_coefficient * normal_iii[2]);

    if (p_i_ii < 0.0)
    {
      if (t1 > 0.0)
      {
        result.region = 4;
        result.stress = {apex, apex, apex};
      }
      else
      {
        result.region = 2;
        result.stress = {t1 + apex, t1 + apex, t1 * friction_coefficient + apex};
      }
    }
    else if (p_i_iii <= 0.0)
    {
      if (t1 > 0.0 && t2 > 0.0)
      {
        result.region = 4;
        result.stress = {apex, apex, apex};
      }
      else
      {
        result.region = 1;
        for (int i = 0; i < 3; ++i)
          result.stress[i] = trial_stress[i] - state_function * return_plane[i];
      }
    }
    else
    {
      if (t2 > 0.0)
      {
        result.region = 4;
        result.stress = {apex, apex, apex};
      }
      else
      {
        result.region = 3;
        result.stress = {
            t2 + apex, t2 * friction_coefficient + apex, t2 * friction_coefficient + apex};
      }
    }

    result.plastic_strain_increment =
        elastic_normal_compliance_times(trial_stress - result.stress, youngs, poisson_ratio);
    return result;
  }

  Core::LinAlg::SymmetricTensor<double, 3, 3> symmetric_square_root(
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& tensor)
  {
    const auto& [eigenvalues, eigenvectors] = Core::LinAlg::eig(tensor);
    Core::LinAlg::SymmetricTensor<double, 3, 3> result{};
    for (int i = 0; i < 3; ++i)
    {
      FOUR_C_ASSERT_ALWAYS(
          eigenvalues[i] > 0.0, "Tensor square root requires positive eigenvalues.");
      Core::LinAlg::Tensor<double, 3> direction{};
      for (int j = 0; j < 3; ++j) direction(j) = eigenvectors(j, i);
      result += std::sqrt(eigenvalues[i]) * Core::LinAlg::self_dyadic(direction);
    }
    return result;
  }
}  // namespace

Mat::PAR::PlasticMohrCoulomb::PlasticMohrCoulomb(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      youngs_(matdata.parameters.get<double>("YOUNG_MODULUS")),
      poisson_ratio_(matdata.parameters.get<double>("POISSON_RATIO")),
      density_(matdata.parameters.get<double>("DENSITY")),
      cohesion_(matdata.parameters.get<double>("COHESION")),
      friction_angle_(matdata.parameters.get<double>("INTERNAL_FRICTION_ANGLE")),
      dilatancy_angle_(matdata.parameters.get<double>("INTERNAL_DILATANCY_ANGLE")),
      isotropic_hardening_(matdata.parameters.get<double>("ISOTROPIC_HARDENING_MODULUS")),
      exponential_saturation_cohesion_(
          matdata.parameters.get<double>("EXPONENTIAL_SATURATION_COHESION")),
      hardening_exponent_(matdata.parameters.get<double>("HARDENING_EXPONENT")),
      tolerance_(matdata.parameters.get<double>("TOLERANCE")),
      max_iterations_(matdata.parameters.get<int>("MAX_ITERATIONS"))
{
  FOUR_C_ASSERT_ALWAYS(exponential_saturation_cohesion_ >= cohesion_,
      "EXPONENTIAL_SATURATION_COHESION must be greater than or equal to COHESION.");
}

std::shared_ptr<Core::Mat::Material> Mat::PAR::PlasticMohrCoulomb::create_material()
{
  return std::make_shared<Mat::PlasticMohrCoulomb>(this);
}

Mat::PlasticMohrCoulombType Mat::PlasticMohrCoulombType::instance_;

Core::Communication::ParObject* Mat::PlasticMohrCoulombType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* material = new Mat::PlasticMohrCoulomb();
  material->unpack(buffer);
  return material;
}

Mat::PlasticMohrCoulomb::PlasticMohrCoulomb() : params_(nullptr) {}

Mat::PlasticMohrCoulomb::PlasticMohrCoulomb(Mat::PAR::PlasticMohrCoulomb* params) : params_(params)
{
}

void Mat::PlasticMohrCoulomb::pack(Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, unique_par_object_id());
  add_to_pack(data, params_ != nullptr ? params_->id() : -1);
  add_to_pack(data, invplrcglast_);
  add_to_pack(data, invplrcgcurr_);
  add_to_pack(data, equivalent_plastic_strain_last_);
  add_to_pack(data, equivalent_plastic_strain_current_);
  add_to_pack(data, accumulated_plastic_volumetric_strain_last_);
  add_to_pack(data, accumulated_plastic_volumetric_strain_current_);
}

void Mat::PlasticMohrCoulomb::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());
  int material_id = -1;
  extract_from_pack(buffer, material_id);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != nullptr &&
      Global::Problem::instance()->materials()->num() != 0)
  {
    const int problem_instance = Global::Problem::instance()->materials()->get_read_from_problem();
    auto* material =
        Global::Problem::instance(problem_instance)->materials()->parameter_by_id(material_id);
    FOUR_C_ASSERT_ALWAYS(material->type() == material_type(),
        "Packed material type does not match MAT_PlasticMohrCoulomb.");
    params_ = static_cast<Mat::PAR::PlasticMohrCoulomb*>(material);
  }
  extract_from_pack(buffer, invplrcglast_);
  extract_from_pack(buffer, invplrcgcurr_);
  extract_from_pack(buffer, equivalent_plastic_strain_last_);
  extract_from_pack(buffer, equivalent_plastic_strain_current_);
  extract_from_pack(buffer, accumulated_plastic_volumetric_strain_last_);
  extract_from_pack(buffer, accumulated_plastic_volumetric_strain_current_);
}

void Mat::PlasticMohrCoulomb::setup(int numgp, const Discret::Elements::Fibers& fibers,
    const std::optional<Discret::Elements::CoordinateSystem>& coord_system)
{
  const auto identity = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;
  invplrcglast_.assign(numgp, identity);
  invplrcgcurr_.assign(numgp, identity);
  equivalent_plastic_strain_last_.assign(numgp, 0.0);
  equivalent_plastic_strain_current_.assign(numgp, 0.0);
  accumulated_plastic_volumetric_strain_last_.assign(numgp, 0.0);
  accumulated_plastic_volumetric_strain_current_.assign(numgp, 0.0);
}

void Mat::PlasticMohrCoulomb::update()
{
  invplrcglast_ = invplrcgcurr_;
  equivalent_plastic_strain_last_ = equivalent_plastic_strain_current_;
  accumulated_plastic_volumetric_strain_last_ = accumulated_plastic_volumetric_strain_current_;
}

double Mat::PlasticMohrCoulomb::cohesion(double equivalent_plastic_strain) const
{
  return params_->cohesion_ + params_->isotropic_hardening_ * equivalent_plastic_strain +
         (params_->exponential_saturation_cohesion_ - params_->cohesion_) *
             (1.0 - std::exp(-params_->hardening_exponent_ * equivalent_plastic_strain));
}

Mat::PlasticMohrCoulomb::ReturnResult Mat::PlasticMohrCoulomb::evaluate_stress(
    const Core::LinAlg::Tensor<double, 3, 3>& deformation_gradient, int gp) const
{
  const auto inverse_deformation_gradient = Core::LinAlg::inv(deformation_gradient);
  const auto elastic_left_cauchy_green_trial = Core::LinAlg::assume_symmetry(
      deformation_gradient * invplrcglast_.at(gp) * Core::LinAlg::transpose(deformation_gradient));
  const auto& [elastic_stretches_squared, eigenvectors] =
      Core::LinAlg::eig(elastic_left_cauchy_green_trial);

  const double shear_modulus = params_->youngs_ / (2.0 * (1.0 + params_->poisson_ratio_));
  const double lame_lambda =
      params_->youngs_ * params_->poisson_ratio_ /
      ((1.0 + params_->poisson_ratio_) * (1.0 - 2.0 * params_->poisson_ratio_));

  Vector3 principal_elastic_strain_unsorted{};
  Vector3 principal_trial_stress_unsorted{};
  for (int i = 0; i < 3; ++i)
  {
    FOUR_C_ASSERT_ALWAYS(elastic_stretches_squared[i] > 0.0,
        "Elastic left Cauchy-Green tensor must be positive definite.");
    principal_elastic_strain_unsorted[i] = 0.5 * std::log(elastic_stretches_squared[i]);
  }
  const double elastic_volume_strain = principal_elastic_strain_unsorted[0] +
                                       principal_elastic_strain_unsorted[1] +
                                       principal_elastic_strain_unsorted[2];
  for (int i = 0; i < 3; ++i)
    principal_trial_stress_unsorted[i] = lame_lambda * elastic_volume_strain +
                                         2.0 * shear_modulus * principal_elastic_strain_unsorted[i];

  std::array<int, 3> order{0, 1, 2};
  std::sort(order.begin(), order.end(), [&](int left, int right)
      { return principal_trial_stress_unsorted[left] > principal_trial_stress_unsorted[right]; });

  Vector3 principal_trial_stress{};
  Vector3 principal_elastic_strain_trial{};
  std::array<Core::LinAlg::Tensor<double, 3>, 3> principal_directions{};
  for (int i = 0; i < 3; ++i)
  {
    principal_trial_stress[i] = principal_trial_stress_unsorted[order[i]];
    principal_elastic_strain_trial[i] = principal_elastic_strain_unsorted[order[i]];
    for (int j = 0; j < 3; ++j) principal_directions[i](j) = eigenvectors(j, order[i]);
  }

  const double old_equivalent_plastic_strain = equivalent_plastic_strain_last_.at(gp);

  auto consistency_residual = [&](double equivalent_plastic_strain)
  {
    const auto returned = return_to_mohr_coulomb(principal_trial_stress,
        cohesion(equivalent_plastic_strain), params_->friction_angle_, params_->dilatancy_angle_,
        params_->youngs_, params_->poisson_ratio_);
    const double equivalent_increment =
        equivalent_plastic_strain_increment(returned.plastic_strain_increment);
    return equivalent_plastic_strain - old_equivalent_plastic_strain - equivalent_increment;
  };

  double new_equivalent_plastic_strain = old_equivalent_plastic_strain;
  if (consistency_residual(old_equivalent_plastic_strain) < 0.0)
  {
    // The stress evaluations used by the numerical tangent must be substantially more accurate
    // than its strain perturbation to avoid force noise in the global Newton iteration.
    const double consistency_tolerance =
        std::min(params_->tolerance_, 1.0e-3 * numerical_tangent_perturbation);
    double lower = old_equivalent_plastic_strain;
    double upper = lower + std::max(1.0e-12, norm(principal_trial_stress) / params_->youngs_);
    int bracket_iteration = 0;
    while (consistency_residual(upper) < 0.0 && bracket_iteration++ < params_->max_iterations_)
      upper = old_equivalent_plastic_strain + 2.0 * (upper - old_equivalent_plastic_strain);
    FOUR_C_ASSERT_ALWAYS(consistency_residual(upper) >= 0.0,
        "Could not bracket the Mohr-Coulomb hardening consistency condition.");

    for (int iteration = 0; iteration < params_->max_iterations_; ++iteration)
    {
      new_equivalent_plastic_strain = 0.5 * (lower + upper);
      const double residual = consistency_residual(new_equivalent_plastic_strain);
      FOUR_C_ASSERT_ALWAYS(std::isfinite(residual),
          "Non-finite Mohr-Coulomb hardening consistency residual at equivalent plastic strain {} "
          "in bracket [{}, {}].",
          new_equivalent_plastic_strain, lower, upper);
      if (std::abs(residual) <=
          consistency_tolerance * std::max(1.0, new_equivalent_plastic_strain))
        break;
      if (residual < 0.0)
        lower = new_equivalent_plastic_strain;
      else
        upper = new_equivalent_plastic_strain;
      FOUR_C_ASSERT_ALWAYS(iteration + 1 < params_->max_iterations_,
          "Mohr-Coulomb hardening consistency solve did not converge in bracket [{}, {}] with "
          "residuals [{}, {}], return regions [{}, {}], and trial principal stress [{}, {}, {}].",
          lower, upper, consistency_residual(lower), consistency_residual(upper),
          return_to_mohr_coulomb(principal_trial_stress, cohesion(lower), params_->friction_angle_,
              params_->dilatancy_angle_, params_->youngs_, params_->poisson_ratio_)
              .region,
          return_to_mohr_coulomb(principal_trial_stress, cohesion(upper), params_->friction_angle_,
              params_->dilatancy_angle_, params_->youngs_, params_->poisson_ratio_)
              .region,
          principal_trial_stress[0], principal_trial_stress[1], principal_trial_stress[2]);
    }
  }

  const auto returned = return_to_mohr_coulomb(principal_trial_stress,
      cohesion(new_equivalent_plastic_strain), params_->friction_angle_, params_->dilatancy_angle_,
      params_->youngs_, params_->poisson_ratio_);
  const Vector3 principal_elastic_strain =
      principal_elastic_strain_trial - returned.plastic_strain_increment;

  ReturnResult result;
  result.equivalent_plastic_strain = new_equivalent_plastic_strain;
  result.accumulated_plastic_volumetric_strain =
      accumulated_plastic_volumetric_strain_last_.at(gp) + returned.plastic_strain_increment[0] +
      returned.plastic_strain_increment[1] + returned.plastic_strain_increment[2];

  for (int i = 0; i < 3; ++i)
  {
    result.elastic_left_cauchy_green += std::exp(2.0 * principal_elastic_strain[i]) *
                                        Core::LinAlg::self_dyadic(principal_directions[i]);
    const auto material_direction = inverse_deformation_gradient * principal_directions[i];
    result.second_piola_kirchhoff_stress +=
        returned.stress[i] * Core::LinAlg::self_dyadic(material_direction);
  }
  result.inverse_plastic_right_cauchy_green = Core::LinAlg::assume_symmetry(
      inverse_deformation_gradient * result.elastic_left_cauchy_green *
      Core::LinAlg::transpose(inverse_deformation_gradient));
  return result;
}

void Mat::PlasticMohrCoulomb::evaluate_numerical_tangent(
    const Core::LinAlg::Tensor<double, 3, 3>& deformation_gradient, int gp,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat) const
{
  const auto right_cauchy_green = Core::LinAlg::assume_symmetry(
      Core::LinAlg::transpose(deformation_gradient) * deformation_gradient);
  const auto right_stretch = symmetric_square_root(right_cauchy_green);
  const auto rotation = deformation_gradient * Core::LinAlg::inv(right_stretch);
  for (int k = 0; k < 3; ++k)
  {
    for (int l = k; l < 3; ++l)
    {
      auto right_cauchy_green_plus = right_cauchy_green;
      auto right_cauchy_green_minus = right_cauchy_green;
      right_cauchy_green_plus(k, l) += 2.0 * numerical_tangent_perturbation;
      right_cauchy_green_minus(k, l) -= 2.0 * numerical_tangent_perturbation;

      const auto deformation_gradient_plus =
          rotation * symmetric_square_root(right_cauchy_green_plus);
      const auto deformation_gradient_minus =
          rotation * symmetric_square_root(right_cauchy_green_minus);
      const auto stress_plus =
          evaluate_stress(deformation_gradient_plus, gp).second_piola_kirchhoff_stress;
      const auto stress_minus =
          evaluate_stress(deformation_gradient_minus, gp).second_piola_kirchhoff_stress;
      const double denominator =
          (k == l) ? 2.0 * numerical_tangent_perturbation : 4.0 * numerical_tangent_perturbation;
      for (int i = 0; i < 3; ++i)
        for (int j = i; j < 3; ++j)
          cmat(i, j, k, l) = (stress_plus(i, j) - stress_minus(i, j)) / denominator;
    }
  }
}

void Mat::PlasticMohrCoulomb::evaluate(
    const Core::LinAlg::Tensor<double, 3, 3>* deformation_gradient,
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
    const Teuchos::ParameterList& parameters, const EvaluationContext<3>& context,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID)
{
  FOUR_C_ASSERT_ALWAYS(
      deformation_gradient != nullptr, "MAT_PlasticMohrCoulomb requires a deformation gradient.");
  const auto result = evaluate_stress(*deformation_gradient, gp);
  stress = result.second_piola_kirchhoff_stress;
  invplrcgcurr_.at(gp) = result.inverse_plastic_right_cauchy_green;
  equivalent_plastic_strain_current_.at(gp) = result.equivalent_plastic_strain;
  accumulated_plastic_volumetric_strain_current_.at(gp) =
      result.accumulated_plastic_volumetric_strain;
  evaluate_numerical_tangent(*deformation_gradient, gp, cmat);
}

Core::LinAlg::SymmetricTensor<double, 3, 3> Mat::PlasticMohrCoulomb::logarithmic_plastic_strain(
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& inverse_plastic_right_cauchy_green) const
{
  const auto& [eigenvalues, eigenvectors] = Core::LinAlg::eig(inverse_plastic_right_cauchy_green);
  Core::LinAlg::SymmetricTensor<double, 3, 3> plastic_strain{};
  for (int i = 0; i < 3; ++i)
  {
    FOUR_C_ASSERT_ALWAYS(eigenvalues[i] > 0.0,
        "Inverse plastic right Cauchy-Green tensor must be positive definite.");
    Core::LinAlg::Tensor<double, 3> direction{};
    for (int j = 0; j < 3; ++j) direction(j) = eigenvectors(j, i);
    plastic_strain += -0.5 * std::log(eigenvalues[i]) * Core::LinAlg::self_dyadic(direction);
  }
  return plastic_strain;
}

void Mat::PlasticMohrCoulomb::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size["plastic_strain"] = 6;
  names_and_size["accumulated_equivalent_plastic_strain"] = 1;
  names_and_size["accumulated_plastic_volumetric_strain"] = 1;
}

bool Mat::PlasticMohrCoulomb::evaluate_output_data(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  if (name == "plastic_strain")
  {
    for (std::size_t gp = 0; gp < invplrcgcurr_.size(); ++gp)
    {
      const auto strain = Core::LinAlg::make_strain_like_voigt_matrix(
          logarithmic_plastic_strain(invplrcgcurr_[gp]));
      for (int component = 0; component < 6; ++component)
        data(gp, component) = strain(component, 0);
    }
    return true;
  }
  if (name == "accumulated_equivalent_plastic_strain")
  {
    for (std::size_t gp = 0; gp < equivalent_plastic_strain_current_.size(); ++gp)
      data(gp, 0) = equivalent_plastic_strain_current_[gp];
    return true;
  }
  if (name == "accumulated_plastic_volumetric_strain")
  {
    for (std::size_t gp = 0; gp < accumulated_plastic_volumetric_strain_current_.size(); ++gp)
      data(gp, 0) = accumulated_plastic_volumetric_strain_current_[gp];
    return true;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
