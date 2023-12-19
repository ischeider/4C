/*---------------------------------------------------------------------------*/
/*! \file
\brief normal contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_H
#define BACI_PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_utils_function_library.H"

#include <Teuchos_ParameterList.hpp>

#include <memory>

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  namespace UTILS
  {
    class FunctionOfScalar;
  }
}  // namespace DRT

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContactNormalBase() = default;

    //! init normal contact handler
    virtual void Init();

    //! setup normal contact handler
    virtual void Setup(const double& dens_max);

    //! get normal contact stiffness
    virtual double GetNormalContactStiffness() const final { return k_normal_; };

    //! get critical normal contact stiffness
    virtual double GetCriticalNormalContactStiffness() const final { return k_normal_crit_; };

    //! get relative penetration
    virtual double GetRelativePenetration() const final { return c_; };

    //! evaluate normal contact force
    virtual void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const = 0;

    //! evaluate normal potential energy
    virtual void NormalPotentialEnergy(const double& gap, double& normalpotentialenergy) const = 0;

   protected:
    //! discrete element method parameter list
    const Teuchos::ParameterList& params_dem_;

    //! maximum expected particle radius
    const double r_max_;

    //! maximum expected particle velocity
    const double v_max_;

    //! maximum allowed relative penetration
    const double c_;

    //! normal contact stiffness
    double k_normal_;

    //! critical normal contact stiffness
    double k_normal_crit_;
  };

  class DEMContactNormalLinearSpring : public virtual DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalLinearSpring(const Teuchos::ParameterList& params);

    //! setup normal contact handler
    void Setup(const double& dens_max) override;

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;

    //! evaluate normal potential energy
    void NormalPotentialEnergy(const double& gap, double& normalpotentialenergy) const override;
  };

  class DEMContactNormalLinearSpringDamp : public DEMContactNormalLinearSpring
  {
   public:
    //! constructor
    explicit DEMContactNormalLinearSpringDamp(const Teuchos::ParameterList& params);

    //! init normal contact handler
    void Init() override;

    //! setup normal contact handler
    void Setup(const double& dens_max) override;

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;

   protected:
    //! coefficient of restitution
    const double e_;

    //! parameter for regularization of damping contact force
    const double damp_reg_fac_;

    //! normal contact damping factor
    double d_normal_fac_;
  };

  class DEMContactNormalNonlinearBase : public DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalNonlinearBase(const Teuchos::ParameterList& params);

    //! setup normal contact handler
    void Setup(const double& dens_max) override;

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override = 0;

    //! evaluate normal potential energy
    void NormalPotentialEnergy(const double& gap, double& normalpotentialenergy) const override;
  };

  class DEMContactNormalHertz : public DEMContactNormalNonlinearBase
  {
   public:
    //! constructor
    explicit DEMContactNormalHertz(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class DEMContactNormalNonlinearDampBase : public DEMContactNormalNonlinearBase
  {
   public:
    //! constructor
    explicit DEMContactNormalNonlinearDampBase(const Teuchos::ParameterList& params);

    //! init normal contact handler
    void Init() override;

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff,
        double& normalcontactforce) const override = 0;

   protected:
    //! normal contact damping parameter
    const double d_normal_;
  };

  class DEMContactNormalLeeHerrmann : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalLeeHerrmann(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class DEMContactNormalKuwabaraKono : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalKuwabaraKono(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class DEMContactNormalTsuji : public DEMContactNormalNonlinearDampBase
  {
   public:
    //! constructor
    explicit DEMContactNormalTsuji(const Teuchos::ParameterList& params);

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;
  };

  class DEMContactNormalForceDistanceCurve : public virtual DEMContactNormalBase
  {
   public:
    //! constructor
    explicit DEMContactNormalForceDistanceCurve(const Teuchos::ParameterList& params);

    //! init normal contact handler
    void Init() override;

    //! setup normal contact handler
    void Setup(const double& dens_max) const;

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;

    //! evaluate normal potential energy
    void NormalPotentialEnergy(const double& gap, double& normalpotentialenergy) const override;

   protected:
    //! curve id which describes force distance relationship
    const int curve_id_;

    //! force distance relationship
    const CORE::UTILS::FunctionOfScalar& force_gap_curve_;
  };

  class DEMContactNormalForceDistanceCurveDamp : public DEMContactNormalForceDistanceCurve,
                                                 public DEMContactNormalLinearSpringDamp
  {
   public:
    //! constructor
    explicit DEMContactNormalForceDistanceCurveDamp(const Teuchos::ParameterList& params);

    //! init normal contact handler
    void Init() override;

    //! setup normal contact handler
    void Setup(const double& dens_max) override;

    //! evaluate normal contact force
    void NormalContactForce(const double& gap, const double radius, const double k_normal,
        const double& v_rel_normal, const double& m_eff, double& normalcontactforce) const override;

    //! evaluate normal potential energy
    void NormalPotentialEnergy(const double& gap, double& normalpotentialenergy) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
