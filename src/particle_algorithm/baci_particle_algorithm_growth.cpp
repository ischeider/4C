/*---------------------------------------------------------------------------*/
/*! \file
\brief particle growth handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_algorithm_growth.H"

#include "baci_lib_globalproblem.H"
#include "baci_particle_engine_container.H"
#include "baci_particle_engine_enums.H"
#include "baci_particle_engine_interface.H"
#include "baci_particle_interaction_utils.H"
#include "baci_utils_function_of_time.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::ParticleGrowthHandler::ParticleGrowthHandler(const int growthfunc)
    : growthfunc_(growthfunc)
{
  // empty constructor
}

void PARTICLEALGORITHM::ParticleGrowthHandler::Init()
{
  // nothing to do
}

void PARTICLEALGORITHM::ParticleGrowthHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEALGORITHM::ParticleGrowthHandler::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // insert states of regular phase particles
    particlestates.insert({PARTICLEENGINE::InitRadius});
  }
}


void PARTICLEALGORITHM::ParticleGrowthHandler::SetInitialStates()
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();
  // iterate over particle types
  for (auto& particleType : particlecontainerbundle->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // scale radius values
    const double* radius = container->GetPtrToState(PARTICLEENGINE::Radius, 0);
    double* initradius = container->GetPtrToState(PARTICLEENGINE::InitRadius, 0);
    for (int i = 0; i < particlestored; ++i)
    {
      initradius[i] = radius[i];
    }
  }
}


void PARTICLEALGORITHM::ParticleGrowthHandler::PreEvaluateTimeStep(double time)
{
  // scale radius and dependent quantities of particles
  if (growthfunc_ > 0)
  {
    const double fac = DRT::Problem::Instance()
                           ->FunctionById<CORE::UTILS::FunctionOfTime>(growthfunc_ - 1)
                           .Evaluate(time);

    if (std::abs(fac - 1.0) < 1e-14) return;

    // get particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
        particleengineinterface_->GetParticleContainerBundle();
    // iterate over particle types
    for (auto& particleType : particlecontainerbundle->GetParticleTypes())
    {
      // get container of owned particles of current particle type
      PARTICLEENGINE::ParticleContainer* container =
          particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

      // get number of particles stored in container
      const int particlestored = container->ParticlesStored();

      // no owned particles of current particle type
      if (particlestored <= 0) continue;

      // scale radius values
      double* radius = container->GetPtrToState(PARTICLEENGINE::Radius, 0);
      double* init_radius = container->GetPtrToState(PARTICLEENGINE::InitRadius, 0);
      double* mass = container->GetPtrToState(PARTICLEENGINE::Mass, 0);
      bool calculate_inertia = container->HaveStoredState(PARTICLEENGINE::Inertia);
      double* inertia;
      if (calculate_inertia) inertia = container->GetPtrToState(PARTICLEENGINE::Inertia, 0);
      for (int i = 0; i < particlestored; ++i)
      {
        double new_radius = init_radius[i] * fac;
        double radius_increase = new_radius / radius[i];
        radius[i] = new_radius;
        // scale mass values
        mass[i] *= PARTICLEINTERACTION::UTILS::Pow<3>(radius_increase);
        // scale inertia values (0.4 mass radius^2)
        if (calculate_inertia) inertia[i] *= PARTICLEINTERACTION::UTILS::Pow<5>(radius_increase);
      }
    }
  }
}


BACI_NAMESPACE_CLOSE
