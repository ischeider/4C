/*---------------------------------------------------------------------------*/
/*! \file
\brief base particle interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_interaction_base.H"

#include "baci_lib_globalproblem.H"
#include "baci_particle_algorithm.H"
#include "baci_particle_engine_container.H"
#include "baci_particle_engine_interface.H"
#include "baci_particle_interaction_material_handler.H"
#include "baci_particle_interaction_runtime_writer.H"
#include "baci_particle_interaction_utils.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionBase::ParticleInteractionBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm), myrank_(comm.MyPID()), params_(params), time_(0.0), dt_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::ParticleInteractionBase::Init()
{
  // init particle material handler
  InitParticleMaterialHandler();

  // init particle interaction writer
  InitParticleInteractionWriter();
}

void PARTICLEINTERACTION::ParticleInteractionBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set interface to particle wall hander
  particlewallinterface_ = particlewallinterface;

  // setup particle material handler
  particlematerial_->Setup();

  // setup particle interaction writer
  particleinteractionwriter_->Setup();

  // init vector
  gravity_.resize(3, 0.0);
}

void PARTICLEINTERACTION::ParticleInteractionBase::WriteRestart() const
{
  // nothing to do
}

void PARTICLEINTERACTION::ParticleInteractionBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of particle interaction writer
  particleinteractionwriter_->ReadRestart(reader);
}

void PARTICLEINTERACTION::ParticleInteractionBase::
    CheckParticleInteractionDistanceConcerningBinSize() const
{
  // get maximum particle interaction distance
  double allprocmaxinteractiondistance = 0.0;
  double maxinteractiondistance = MaxInteractionDistance();
  comm_.MaxAll(&maxinteractiondistance, &allprocmaxinteractiondistance, 1);

  // bin size safety check
  if (allprocmaxinteractiondistance > particleengineinterface_->MinBinSize())
    dserror("the particle interaction distance is larger than the minimal bin size (%f > %f)!",
        allprocmaxinteractiondistance, particleengineinterface_->MinBinSize());

  // periodic length safety check
  if (particleengineinterface_->HavePeriodicBoundaryConditions())
  {
    // loop over all spatial directions
    for (int dim = 0; dim < 3; ++dim)
    {
      // check for periodic boundary condition in current spatial direction
      if (not particleengineinterface_->HavePeriodicBoundaryConditionsInSpatialDirection(dim))
        continue;

      // check periodic length in current spatial direction
      if ((2.0 * allprocmaxinteractiondistance) >
          particleengineinterface_->LengthOfBinningDomainInASpatialDirection(dim))
        dserror("particles are not allowed to interact directly and across the periodic boundary!");
    }
  }
}

void PARTICLEINTERACTION::ParticleInteractionBase::SetCurrentTime(const double currenttime)
{
  time_ = currenttime;
}

void PARTICLEINTERACTION::ParticleInteractionBase::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

void PARTICLEINTERACTION::ParticleInteractionBase::SetCurrentWriteResultFlag(
    bool writeresultsthisstep)
{
  // set current write result flag in particle interaction writer
  particleinteractionwriter_->SetCurrentWriteResultFlag(writeresultsthisstep);
}

void PARTICLEINTERACTION::ParticleInteractionBase::SetGravity(std::vector<double>& gravity)
{
  gravity_ = gravity;
}

void PARTICLEINTERACTION::ParticleInteractionBase::WriteInteractionRuntimeOutput(
    const int step, const double time)
{
  // write particle interaction runtime output
  particleinteractionwriter_->WriteParticleInteractionRuntimeOutput(step, time);
}

void PARTICLEINTERACTION::ParticleInteractionBase::InitParticleMaterialHandler()
{
  // create particle material handler
  particlematerial_ = std::make_shared<PARTICLEINTERACTION::MaterialHandler>(params_);

  // init particle material handler
  particlematerial_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionBase::InitParticleInteractionWriter()
{
  // create particle interaction writer
  particleinteractionwriter_ =
      std::make_shared<PARTICLEINTERACTION::InteractionWriter>(comm_, params_);

  // init particle interaction writer
  particleinteractionwriter_->Init();
}

double PARTICLEINTERACTION::ParticleInteractionBase::MaxParticleRadius() const
{
  // init value of maximum radius
  double maxrad = 0.0;

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get maximum stored value of state
    double currmaxrad = container->GetMaxValueOfState(PARTICLEENGINE::Radius);

    // compare to current maximum
    maxrad = std::max(maxrad, currmaxrad);
  }

  return maxrad;
}

BACI_NAMESPACE_CLOSE
