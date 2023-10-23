/*---------------------------------------------------------------------------*/
/*! \file
\brief discrete element method (DEM) interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_DEM_H
#define BACI_PARTICLE_INTERACTION_DEM_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_inpar_particle.H"
#include "baci_particle_interaction_base.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMNeighborPairs;
  class DEMHistoryPairs;
  class DEMContact;
  class DEMAdhesion;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  /*!
   * \brief discrete element method (DEM) interaction
   *
   * \author Sebastian Fuchs \date 05/2018
   */

  class ParticleInteractionDEM final : public ParticleInteractionBase
  {
   public:
    //! constructor
    explicit ParticleInteractionDEM(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 10/2018
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~ParticleInteractionDEM() override;

    //! init particle interaction handler
    void Init() override;

    //! setup particle interaction handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface) override;

    //! write restart of particle interaction handler
    void WriteRestart() const override;

    //! read restart of particle interaction handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader) override;

    //! insert interaction dependent states of all particle types
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) override;

    //! set initial states
    void SetInitialStates() override;

    //! pre evaluate time step
    void PreEvaluateTimeStep() override;

    //! evaluate particle interactions
    void EvaluateInteractions() override;

    //! post evaluate time step
    void PostEvaluateTimeStep(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) override;

    //! maximum interaction distance (on this processor)
    double MaxInteractionDistance() const override;

    //! distribute interaction history
    void DistributeInteractionHistory() const override;

    //! communicate interaction history
    void CommunicateInteractionHistory() const override;

    //! set current step size
    void SetCurrentStepSize(const double currentstepsize) override;

   private:
    //! init neighbor pair handler
    void InitNeighborPairHandler();

    //! init history pair handler
    void InitHistoryPairHandler();

    //! init contact handler
    void InitContactHandler();

    //! init adhesion handler
    void InitAdhesionHandler();

    //! setup particle interaction writer
    void SetupParticleInteractionWriter();

    //! set initial radius
    void SetInitialRadius();

    //! set initial mass
    void SetInitialMass();

    //! set initial inertia
    void SetInitialInertia();

    //! clear force and moment states of particles
    void ClearForceAndMomentStates() const;

    //! compute acceleration from force and moment
    void ComputeAcceleration() const;

    //! evaluate particle energy
    void EvaluateParticleEnergy() const;

    //! evaluate particle kinetic energy contribution
    void EvaluateParticleKineticEnergy(double& kineticenergy) const;

    //! evaluate particle gravitational potential energy contribution
    void EvaluateParticleGravitationalPotentialEnergy(double& gravitationalpotentialenergy) const;

    //! discrete element method specific parameter list
    const Teuchos::ParameterList& params_dem_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs_;

    //! history pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs_;

    //! contact handler
    std::unique_ptr<PARTICLEINTERACTION::DEMContact> contact_;

    //! adhesion handler
    std::unique_ptr<PARTICLEINTERACTION::DEMAdhesion> adhesion_;

    //! write particle energy output
    const bool writeparticleenergy_;

    //! if gap is smaller than this distance, interaction forces are evaluated
    const double interaction_horizon_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
