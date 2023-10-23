/*---------------------------------------------------------------------------*/
/*! \file
\brief normal contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_interaction_dem_contact_normal.H"

#include "baci_inpar_particle.H"
#include "baci_lib_function_library.H"
#include "baci_lib_globalproblem.H"
#include "baci_particle_interaction_utils.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContactNormalBase::DEMContactNormalBase(
    const Teuchos::ParameterList& params)
    : params_dem_(params),
      r_max_(params_dem_.get<double>("MAX_RADIUS")),
      v_max_(params_dem_.get<double>("MAX_VELOCITY")),
      c_(params_dem_.get<double>("REL_PENETRATION")),
      k_normal_(params_dem_.get<double>("NORMAL_STIFF")),
      k_normal_crit_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalBase::Init()
{
  if (not((c_ <= 0.0 and k_normal_ > 0.0) or (c_ > 0.0 and v_max_ > 0.0 and k_normal_ <= 0.0)))
    dserror(
        "specify either the relative penetration along with the maximum velocity, or the normal "
        "stiffness, but neither both nor none of them!");
}

void PARTICLEINTERACTION::DEMContactNormalBase::Setup(const double& dens_max)
{
  // nothing to do
}

PARTICLEINTERACTION::DEMContactNormalLinearSpring::DEMContactNormalLinearSpring(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalLinearSpring::Setup(const double& dens_max)
{
  // call base class setup
  DEMContactNormalBase::Setup(dens_max);

  // calculate normal stiffness from relative penetration and other input parameters
  if (c_ > 0.0)
    k_normal_ = 2.0 / 3.0 * r_max_ * M_PI * dens_max * UTILS::Pow<2>(v_max_) / UTILS::Pow<2>(c_);

  // set critical normal contact stiffness to linear normal contact stiffness
  k_normal_crit_ = k_normal_;
}

void PARTICLEINTERACTION::DEMContactNormalLinearSpring::NormalContactForce(const double& gap,
    const double* radius_i, const double* radius_j, const double& v_rel_normal, const double& m_eff,
    double& normalcontactforce) const
{
  normalcontactforce = k_normal_ * gap;
}

void PARTICLEINTERACTION::DEMContactNormalLinearSpring::NormalPotentialEnergy(
    const double& gap, double& normalpotentialenergy) const
{
  normalpotentialenergy = 0.5 * k_normal_ * UTILS::Pow<2>(gap);
}

PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp::DEMContactNormalLinearSpringDamp(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalBase(params),
      PARTICLEINTERACTION::DEMContactNormalLinearSpring(params),
      e_(params_dem_.get<double>("COEFF_RESTITUTION")),
      damp_reg_fac_(params_dem_.get<double>("DAMP_REG_FAC")),
      d_normal_fac_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp::Init()
{
  // call base class init
  DEMContactNormalLinearSpring::Init();

  // safety checks for contact parameters
  if (e_ < 0.0) dserror("invalid input parameter COEFF_RESTITUTION for this kind of contact law!");
}

void PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp::Setup(const double& dens_max)
{
  // call base class setup
  DEMContactNormalLinearSpring::Setup(dens_max);

  // determine normal contact damping factor
  if (e_ > 0.0)
  {
    const double lne = std::log(e_);
    d_normal_fac_ =
        2.0 * std::abs(lne) * std::sqrt(k_normal_ / (UTILS::Pow<2>(lne) + UTILS::Pow<2>(M_PI)));
  }
  else
    d_normal_fac_ = 2.0 * std::sqrt(k_normal_);
}

void PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp::NormalContactForce(const double& gap,
    const double* radius_i, const double* radius_j, const double& v_rel_normal, const double& m_eff,
    double& normalcontactforce) const
{
  // determine normal contact damping parameter
  double d_normal = d_normal_fac_ * std::sqrt(m_eff);

  // linear regularization of damping force to reach full amplitude for g = damp_reg_fac_ * r_min
  double reg_fac = 1.0;
  if (damp_reg_fac_ > 0.0)
  {
    double rad_min = std::min(radius_i[0], radius_j[0]);

    if (std::abs(gap) < damp_reg_fac_ * rad_min)
      reg_fac = (std::abs(gap) / (damp_reg_fac_ * rad_min));
  }

  // evaluate normal contact force
  normalcontactforce = (k_normal_ * gap) - (d_normal * v_rel_normal * reg_fac);
}

PARTICLEINTERACTION::DEMContactNormalNonlinearBase::DEMContactNormalNonlinearBase(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalNonlinearBase::Setup(const double& dens_max)
{
  // call base class setup
  DEMContactNormalBase::Setup(dens_max);

  // calculate normal stiffness from relative penetration and other input parameters if necessary
  if (c_ > 0.0)
    k_normal_ = 10.0 / 3.0 * M_PI * dens_max * UTILS::Pow<2>(v_max_) * std::sqrt(r_max_) /
                std::sqrt(UTILS::Pow<5>(2.0 * c_));

  // set critical normal contact stiffness to linear normal contact stiffness
  if (c_ > 0.0)
    k_normal_crit_ =
        2.0 / 3.0 * r_max_ * M_PI * dens_max * UTILS::Pow<2>(v_max_) / UTILS::Pow<2>(c_);
  else
    k_normal_crit_ = std::pow(2048.0 / 1875.0 * dens_max * UTILS::Pow<2>(v_max_) * M_PI *
                                  UTILS::Pow<3>(r_max_) * UTILS::Pow<4>(k_normal_),
        0.2);
}

void PARTICLEINTERACTION::DEMContactNormalNonlinearBase::NormalPotentialEnergy(
    const double& gap, double& normalpotentialenergy) const
{
  normalpotentialenergy = 0.4 * k_normal_ * UTILS::Pow<2>(gap) * std::sqrt(-gap);
}

PARTICLEINTERACTION::DEMContactNormalHertz::DEMContactNormalHertz(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalNonlinearBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalHertz::NormalContactForce(const double& gap,
    const double* radius_i, const double* radius_j, const double& v_rel_normal, const double& m_eff,
    double& normalcontactforce) const
{
  normalcontactforce = -k_normal_ * (-gap) * std::sqrt(-gap);
}

PARTICLEINTERACTION::DEMContactNormalNonlinearDampBase::DEMContactNormalNonlinearDampBase(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalNonlinearBase(params),
      d_normal_(params_dem_.get<double>("NORMAL_DAMP"))
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalNonlinearDampBase::Init()
{
  // call base class init
  DEMContactNormalNonlinearBase::Init();

  // safety checks for contact parameters
  if (d_normal_ < 0.0) dserror("invalid input parameter NORMAL_DAMP for this kind of contact law!");
}

PARTICLEINTERACTION::DEMContactNormalLeeHerrmann::DEMContactNormalLeeHerrmann(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalNonlinearDampBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalLeeHerrmann::NormalContactForce(const double& gap,
    const double* radius_i, const double* radius_j, const double& v_rel_normal, const double& m_eff,
    double& normalcontactforce) const
{
  normalcontactforce = -k_normal_ * (-gap) * std::sqrt(-gap) - m_eff * d_normal_ * v_rel_normal;
}

PARTICLEINTERACTION::DEMContactNormalKuwabaraKono::DEMContactNormalKuwabaraKono(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalNonlinearDampBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalKuwabaraKono::NormalContactForce(const double& gap,
    const double* radius_i, const double* radius_j, const double& v_rel_normal, const double& m_eff,
    double& normalcontactforce) const
{
  normalcontactforce =
      -k_normal_ * (-gap) * std::sqrt(-gap) - d_normal_ * v_rel_normal * std::sqrt(-gap);
}

PARTICLEINTERACTION::DEMContactNormalTsuji::DEMContactNormalTsuji(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalNonlinearDampBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalTsuji::NormalContactForce(const double& gap,
    const double* radius_i, const double* radius_j, const double& v_rel_normal, const double& m_eff,
    double& normalcontactforce) const
{
  normalcontactforce =
      -k_normal_ * (-gap) * std::sqrt(-gap) - d_normal_ * v_rel_normal * std::sqrt(std::sqrt(-gap));
}

PARTICLEINTERACTION::DEMContactNormalForceDistanceCurve::DEMContactNormalForceDistanceCurve(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalBase(params),
      curve_id_(params_dem_.get<int>("NORMAL_CONTACT_FORCE_DIST_CURVE")),
      force_gap_curve_(
          DRT::Problem::Instance()->FunctionById<DRT::UTILS::FunctionOfScalar>(curve_id_ - 1))
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalForceDistanceCurve::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMContactNormalForceDistanceCurve::Setup(const double& dens_max) const
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMContactNormalForceDistanceCurve::NormalContactForce(const double& gap,
    const double* radius_i, const double* radius_j, const double& v_rel_normal, const double& m_eff,
    double& normalcontactforce) const
{
  normalcontactforce = force_gap_curve_.Evaluate(gap);
}

void PARTICLEINTERACTION::DEMContactNormalForceDistanceCurve::NormalPotentialEnergy(
    const double& gap, double& normalpotentialenergy) const
{
  normalpotentialenergy = 0.5 * force_gap_curve_.EvaluateDerivative(gap, 1) * UTILS::Pow<2>(gap);
}

PARTICLEINTERACTION::DEMContactNormalForceDistanceCurveDamp::DEMContactNormalForceDistanceCurveDamp(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMContactNormalBase(params),
      PARTICLEINTERACTION::DEMContactNormalForceDistanceCurve(params),
      PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMContactNormalForceDistanceCurveDamp::Init()
{
  // call base class init
  DEMContactNormalLinearSpringDamp::Init();
}

void PARTICLEINTERACTION::DEMContactNormalForceDistanceCurveDamp::Setup(const double& dens_max)
{
  // call base class setup
  DEMContactNormalLinearSpringDamp::Setup(dens_max);
}

void PARTICLEINTERACTION::DEMContactNormalForceDistanceCurveDamp::NormalContactForce(
    const double& gap, const double* radius_i, const double* radius_j, const double& v_rel_normal,
    const double& m_eff, double& normalcontactforce) const
{
  // determine normal contact damping parameter
  double d_normal = d_normal_fac_ * std::sqrt(m_eff);

  // linear regularization of damping force to reach full amplitude for g = damp_reg_fac_ * r_min
  double reg_fac = 1.0;
  if (damp_reg_fac_ > 0.0)
  {
    double rad_min = std::min(radius_i[0], radius_j[0]);

    if (std::abs(gap) < damp_reg_fac_ * rad_min)
      reg_fac = (std::abs(gap) / (damp_reg_fac_ * rad_min));
  }

  // evaluate normal contact force
  double d1 = force_gap_curve_.Evaluate(gap);
  double d2 = (d_normal * v_rel_normal * reg_fac);
  normalcontactforce = force_gap_curve_.Evaluate(gap) - (d_normal * v_rel_normal * reg_fac);
}

void PARTICLEINTERACTION::DEMContactNormalForceDistanceCurveDamp::NormalPotentialEnergy(
    const double& gap, double& normalpotentialenergy) const
{
  // call base class energy evaluation
  DEMContactNormalForceDistanceCurve::NormalPotentialEnergy(gap, normalpotentialenergy);
}

BACI_NAMESPACE_CLOSE
