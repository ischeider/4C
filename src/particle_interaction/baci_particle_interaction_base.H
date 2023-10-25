/*---------------------------------------------------------------------------*/
/*! \file
\brief base particle interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_BASE_H
#define BACI_PARTICLE_INTERACTION_BASE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_linalg_fixedsizematrix.H"
#include "baci_particle_engine_enums.H"
#include "baci_particle_engine_typedefs.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEINTERACTION
{
  class MaterialHandler;
  class InteractionWriter;
}  // namespace PARTICLEINTERACTION

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

namespace PARTICLEWALL
{
  class WallHandlerInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  /*!
   * \brief base particle interaction
   *
   * \author Sebastian Fuchs \date 05/2018
   */
  class ParticleInteractionBase
  {
   public:
    //! constructor
    explicit ParticleInteractionBase(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~ParticleInteractionBase() = default;

    //! init particle interaction handler
    virtual void Init();

    //! setup particle interaction handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface);

    //! write restart of particle interaction handler
    virtual void WriteRestart() const;

    //! read restart of particle interaction handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! insert interaction dependent states of all particle types
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) = 0;

    //! set initial states
    virtual void SetInitialStates() = 0;

    //! pre evaluate time step
    virtual void PreEvaluateTimeStep() = 0;

    //! evaluate particle interactions
    virtual void EvaluateInteractions() = 0;

    //! post evaluate time step
    virtual void PostEvaluateTimeStep(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) = 0;

    //! check particle interaction distance concerning bin size
    virtual void CheckParticleInteractionDistanceConcerningBinSize() const final;

    //! maximum interaction distance (on this processor)
    virtual double MaxInteractionDistance() const = 0;

    //! distribute interaction history
    virtual void DistributeInteractionHistory() const = 0;

    //! communicate interaction history
    virtual void CommunicateInteractionHistory() const = 0;

    //! set current time
    virtual void SetCurrentTime(const double currenttime);

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize);

    //! set current write result flag
    virtual void SetCurrentWriteResultFlag(bool writeresultsthisstep);

    //! set gravity
    virtual void SetGravity(std::vector<double>& gravity) final;

    //! write interaction runtime output
    virtual void WriteInteractionRuntimeOutput(const int step, const double time);

   private:
    //! init particle material handler
    virtual void InitParticleMaterialHandler();

    //! init particle interaction writer
    virtual void InitParticleInteractionWriter();

   protected:
    //! maximum particle radius (on this processor)
    virtual double MaxParticleRadius() const;

    //! communication
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

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

    //! current time
    double time_;

    //! time step size
    double dt_;

    //! current gravity vector
    std::vector<double> gravity_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
