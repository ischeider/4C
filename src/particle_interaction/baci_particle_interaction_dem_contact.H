/*---------------------------------------------------------------------------*/
/*! \file
\brief contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_DEM_CONTACT_H
#define BACI_PARTICLE_INTERACTION_DEM_CONTACT_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_inpar_particle.H"
#include "baci_particle_engine_enums.H"
#include "baci_particle_engine_typedefs.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEWALL
{
  class WallHandlerInterface;
}

namespace PARTICLEINTERACTION
{
  class MaterialHandler;
  class InteractionWriter;
  class DEMNeighborPairs;
  class DEMHistoryPairs;
  class DEMContactNormalBase;
  class DEMContactTangentialBase;
  class DEMContactRollingBase;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContact final
  {
   public:
    //! constructor
    explicit DEMContact(const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~DEMContact();

    //! init contact handler
    void Init();

    //! setup contact handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter,
        const std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs,
        const std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs);

    //! set current step size
    void SetCurrentStepSize(const double currentstepsize);

    //! insert contact evaluation dependent states
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! get normal contact stiffness
    double GetNormalContactStiffness() const;

    //! check critical time step (on this processor)
    double CheckCriticalTimeStep() const;

    //! add contact contribution to force and moment field
    void AddForceAndMomentContribution();

    //! evaluate elastic potential energy contribution
    void EvaluateElasticPotentialEnergy(double& elasticpotentialenergy) const;

   private:
    //! init normal contact handler
    void InitNormalContactHandler();

    //! init tangential contact handler
    void InitTangentialContactHandler();

    //! init rolling contact handler
    void InitRollingContactHandler();

    //! setup particle interaction writer
    void SetupParticleInteractionWriter();

    //! get maximum density of all materials
    double GetMaxDensityOfAllMaterials() const;

    //! evaluate particle contact contribution
    void EvaluateParticleContact();

    //! evaluate particle-wall contact contribution
    void EvaluateParticleWallContact();

    //! evaluate particle elastic potential energy contribution
    void EvaluateParticleElasticPotentialEnergy(double& elasticpotentialenergy) const;

    //! evaluate particle-wall elastic potential energy contribution
    void EvaluateParticleWallElasticPotentialEnergy(double& elasticpotentialenergy) const;

    //! discrete element method specific parameter list
    const Teuchos::ParameterList& params_dem_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! particle interaction writer
    std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs_;

    //! history pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs_;

    //! normal contact handler
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalBase> contactnormal_;

    //! tangential contact handler
    std::unique_ptr<PARTICLEINTERACTION::DEMContactTangentialBase> contacttangential_;

    //! rolling contact handler
    std::unique_ptr<PARTICLEINTERACTION::DEMContactRollingBase> contactrolling_;

    //! time step size
    double dt_;

    //! tension cutoff of normal contact force
    const bool tension_cutoff_;

    //! write particle-wall interaction output
    const bool writeparticlewallinteraction_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
