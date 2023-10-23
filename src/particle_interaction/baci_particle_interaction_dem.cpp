/*---------------------------------------------------------------------------*/
/*! \file
\brief discrete element method (DEM) interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_interaction_dem.H"

#include "baci_io_runtime_csv_writer.H"
#include "baci_lib_globalproblem.H"
#include "baci_particle_engine_container.H"
#include "baci_particle_engine_interface.H"
#include "baci_particle_interaction_dem_adhesion.H"
#include "baci_particle_interaction_dem_contact.H"
#include "baci_particle_interaction_dem_history_pairs.H"
#include "baci_particle_interaction_dem_neighbor_pairs.H"
#include "baci_particle_interaction_material_handler.H"
#include "baci_particle_interaction_runtime_writer.H"
#include "baci_particle_interaction_utils.H"
#include "baci_particle_wall_interface.H"

#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionDEM::ParticleInteractionDEM(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::ParticleInteractionBase(comm, params),
      params_dem_(params.sublist("DEM")),
      writeparticleenergy_(DRT::INPUT::IntegralValue<int>(params_dem_, "WRITE_PARTICLE_ENERGY")),
      interaction_horizon_(params_dem_.get<double>("NORMAL_CONTACT_INTERACTION_HORIZON"))
{
  // empty constructor
}

PARTICLEINTERACTION::ParticleInteractionDEM::~ParticleInteractionDEM() = default;

void PARTICLEINTERACTION::ParticleInteractionDEM::Init()
{
  // call base class init
  ParticleInteractionBase::Init();

  // init neighbor pair handler
  InitNeighborPairHandler();

  // init history pair handler
  InitHistoryPairHandler();

  // init contact handler
  InitContactHandler();

  // init adhesion handler
  InitAdhesionHandler();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // call base class setup
  ParticleInteractionBase::Setup(particleengineinterface, particlewallinterface);

  // setup neighbor pair handler
  neighborpairs_->Setup(particleengineinterface, particlewallinterface);

  // setup history pair handler
  historypairs_->Setup(particleengineinterface);

  // setup contact handler
  contact_->Setup(particleengineinterface, particlewallinterface, particlematerial_,
      particleinteractionwriter_, neighborpairs_, historypairs_);

  // setup adhesion handler
  if (adhesion_)
    adhesion_->Setup(particleengineinterface, particlewallinterface, particleinteractionwriter_,
        neighborpairs_, historypairs_, contact_->GetNormalContactStiffness());

  // setup particle interaction writer
  SetupParticleInteractionWriter();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::WriteRestart() const
{
  // call base class function
  ParticleInteractionBase::WriteRestart();

  // write restart of history pair handler
  historypairs_->WriteRestart();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // call base class function
  ParticleInteractionBase::ReadRestart(reader);

  // read restart of history pair handler
  historypairs_->ReadRestart(reader);
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // insert states of regular phase particles
    particlestates.insert({PARTICLEENGINE::Force, PARTICLEENGINE::Mass, PARTICLEENGINE::Radius});
  }

  // states for contact evaluation scheme
  contact_->InsertParticleStatesOfParticleTypes(particlestatestotypes);
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialStates()
{
  // set initial radius
  SetInitialRadius();

  // set initial mass
  SetInitialMass();

  // set initial inertia
  SetInitialInertia();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::PreEvaluateTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::PreEvaluateTimeStep");

  // call base class pre evaluate
  ParticleInteractionBase::PreEvaluateTimeStep();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateInteractions()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateInteractions");

  // clear force and moment states of particles
  ClearForceAndMomentStates();

  // evaluate neighbor pairs
  neighborpairs_->EvaluateNeighborPairs(interaction_horizon_);

  // evaluate adhesion neighbor pairs
  if (adhesion_) neighborpairs_->EvaluateNeighborPairsAdhesion(adhesion_->GetAdhesionDistance());

  // check critical time step
  contact_->CheckCriticalTimeStep();

  // add contact contribution to force and moment field
  contact_->AddForceAndMomentContribution();

  // add adhesion contribution to force field
  if (adhesion_) adhesion_->AddForceContribution();

  // compute acceleration from force and moment
  ComputeAcceleration();

  // update history pairs
  historypairs_->UpdateHistoryPairs();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::PostEvaluateTimeStep(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::PostEvaluateTimeStep");

  // evaluate particle energy
  if (particleinteractionwriter_->GetCurrentWriteResultFlag() and writeparticleenergy_)
    EvaluateParticleEnergy();
}

double PARTICLEINTERACTION::ParticleInteractionDEM::MaxInteractionDistance() const
{
  // particle contact interaction distance
  double interactiondistance = 2.0 * MaxParticleRadius();

  // add adhesion distance
  if (adhesion_) interactiondistance += adhesion_->GetAdhesionDistance();

  return interactiondistance;
}

void PARTICLEINTERACTION::ParticleInteractionDEM::DistributeInteractionHistory() const
{
  // distribute history pairs
  historypairs_->DistributeHistoryPairs();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::CommunicateInteractionHistory() const
{
  // communicate history pairs
  historypairs_->CommunicateHistoryPairs();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetCurrentStepSize(const double currentstepsize)
{
  // call base class method
  ParticleInteractionBase::SetCurrentStepSize(currentstepsize);

  // set current step size
  contact_->SetCurrentStepSize(currentstepsize);
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InitNeighborPairHandler()
{
  // create neighbor pair handler
  neighborpairs_ = std::make_shared<PARTICLEINTERACTION::DEMNeighborPairs>();

  // init neighbor pair handler
  neighborpairs_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InitHistoryPairHandler()
{
  // create history pair handler
  historypairs_ = std::make_shared<PARTICLEINTERACTION::DEMHistoryPairs>(comm_);

  // init history pair handler
  historypairs_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InitContactHandler()
{
  // create contact handler
  contact_ = std::unique_ptr<PARTICLEINTERACTION::DEMContact>(
      new PARTICLEINTERACTION::DEMContact(params_dem_));

  // init contact handler
  contact_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::InitAdhesionHandler()
{
  // get type of adhesion law
  INPAR::PARTICLE::AdhesionLaw adhesionlaw =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::AdhesionLaw>(params_dem_, "ADHESIONLAW");

  // create adhesion handler
  if (adhesionlaw != INPAR::PARTICLE::NoAdhesion)
    adhesion_ = std::unique_ptr<PARTICLEINTERACTION::DEMAdhesion>(
        new PARTICLEINTERACTION::DEMAdhesion(params_dem_));

  // init adhesion handler
  if (adhesion_) adhesion_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetupParticleInteractionWriter()
{
  if (writeparticleenergy_)
  {
    // register specific runtime csv writer
    particleinteractionwriter_->RegisterSpecificRuntimeCsvWriter("particle-energy");

    // get specific runtime csv writer
    IO::RuntimeCsvWriter* runtime_csv_writer =
        particleinteractionwriter_->GetSpecificRuntimeCsvWriter("particle-energy");

    // register all data vectors
    runtime_csv_writer->RegisterDataVector("kin_energy", 1, 10);
    runtime_csv_writer->RegisterDataVector("grav_pot_energy", 1, 10);
    runtime_csv_writer->RegisterDataVector("elast_pot_energy", 1, 10);
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialRadius()
{
  // get allowed bounds for particle radius
  double r_min = params_dem_.get<double>("MIN_RADIUS");
  double r_max = params_dem_.get<double>("MAX_RADIUS");

  // safety checks
  if (r_min < 0.0) dserror("negative minimum allowed particle radius!");
  if (not(r_max > 0.0)) dserror("non-positive maximum allowed particle radius!");
  if (r_min > r_max)
    dserror("minimum allowed particle radius larger than maximum allowed particle radius!");

  // get type of initial particle radius assignment
  INPAR::PARTICLE::InitialRadiusAssignment radiusdistributiontype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::InitialRadiusAssignment>(
          params_dem_, "INITIAL_RADIUS");

  switch (radiusdistributiontype)
  {
    // particle radius from particle material
    case INPAR::PARTICLE::RadiusFromParticleMaterial:
    {
      // iterate over particle types
      for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
      {
        // get container of owned particles of current particle type
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

        // get number of particles stored in container
        const int particlestored = container->ParticlesStored();

        // no owned particles of current particle type
        if (particlestored <= 0) continue;

        // get material for current particle type
        const MAT::PAR::ParticleMaterialBase* material =
            particlematerial_->GetPtrToParticleMatParameter(type_i);

        // safety checks
        if (material->initRadius_ < r_min)
          dserror("material particle radius smaller than minimum allowed particle radius!");

        if (material->initRadius_ > r_max)
          dserror("material particle radius larger than maximum allowed particle radius!");

        // (initial) radius of current phase
        std::vector<double> initradius(1);
        initradius[0] = material->initRadius_;

        // set initial radius for all particles of current type
        container->SetState(initradius, PARTICLEENGINE::Radius);
      }

      break;
    }
    // particle radius from particle input
    case INPAR::PARTICLE::RadiusFromParticleInput:
    {
      // note: particle radius set as read in from input file, only safety checks here

      // iterate over particle types
      for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
      {
        // get container of owned particles of current particle type
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

        // get number of particles stored in container
        const int particlestored = container->ParticlesStored();

        // no owned particles of current particle type
        if (particlestored <= 0) continue;

        // safety checks
        if (container->GetMinValueOfState(PARTICLEENGINE::Radius) < r_min)
          dserror("minimum particle radius smaller than minimum allowed particle radius!");

        if (container->GetMaxValueOfState(PARTICLEENGINE::Radius) > r_max)
          dserror("maximum particle radius larger than maximum allowed particle radius!");
      }

      break;
    }
    // normal or log-normal random particle radius distribution
    case INPAR::PARTICLE::NormalRadiusDistribution:
    case INPAR::PARTICLE::LogNormalRadiusDistribution:
    {
      // get sigma of random particle radius distribution
      double sigma = params_dem_.get<double>("RADIUSDISTRIBUTION_SIGMA");

      // safety check
      if (not(sigma > 0.0)) dserror("non-positive sigma of random particle radius distribution!");

      // iterate over particle types
      for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
      {
        // get container of owned particles of current particle type
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

        // get number of particles stored in container
        const int particlestored = container->ParticlesStored();

        // no owned particles of current particle type
        if (particlestored <= 0) continue;

        // get material for current particle type
        const MAT::PAR::ParticleMaterialBase* material =
            particlematerial_->GetPtrToParticleMatParameter(type_i);

        // get pointer to particle state
        double* radius = container->GetPtrToState(PARTICLEENGINE::Radius, 0);

        // determine mu of random particle radius distribution
        const double mu = (radiusdistributiontype == INPAR::PARTICLE::NormalRadiusDistribution)
                              ? material->initRadius_
                              : std::log(material->initRadius_);

        // initialize random number generator
        DRT::Problem::Instance()->Random()->SetMeanVariance(mu, sigma);

        // iterate over particles stored in container
        for (int i = 0; i < particlestored; ++i)
        {
          // generate random value
          const double randomvalue = DRT::Problem::Instance()->Random()->Normal();

          // set normal or log-normal distributed random value for particle radius
          radius[i] = (radiusdistributiontype == INPAR::PARTICLE::NormalRadiusDistribution)
                          ? randomvalue
                          : std::exp(randomvalue);

          // adjust radius to allowed bounds
          if (radius[i] > r_max)
            radius[i] = r_max;
          else if (radius[i] < r_min)
            radius[i] = r_min;
        }
      }
      break;
    }
    default:
    {
      dserror("invalid type of (random) particle radius distribution!");
      break;
    }
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialMass()
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // get pointer to particle states
    const double* radius = container->GetPtrToState(PARTICLEENGINE::Radius, 0);
    double* mass = container->GetPtrToState(PARTICLEENGINE::Mass, 0);

    // compute mass via particle volume and initial density
    const double fac = material->initDensity_ * 4.0 / 3.0 * M_PI;
    for (int i = 0; i < particlestored; ++i) mass[i] = fac * UTILS::Pow<3>(radius[i]);
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::SetInitialInertia()
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // no inertia state for current particle type
    if (not container->HaveStoredState(PARTICLEENGINE::Inertia)) continue;

    // get pointer to particle states
    const double* radius = container->GetPtrToState(PARTICLEENGINE::Radius, 0);
    const double* mass = container->GetPtrToState(PARTICLEENGINE::Mass, 0);
    double* inertia = container->GetPtrToState(PARTICLEENGINE::Inertia, 0);

    // inertia = 2/5 mass times radius squared
    for (int i = 0; i < particlestored; ++i) inertia[i] = 0.4 * mass[i] * UTILS::Pow<2>(radius[i]);
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::ClearForceAndMomentStates() const
{
  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear force of all particles
    container->ClearState(PARTICLEENGINE::Force);

    // clear moment of all particles
    if (container->HaveStoredState(PARTICLEENGINE::Moment))
      container->ClearState(PARTICLEENGINE::Moment);
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::ComputeAcceleration() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::ComputeAcceleration");

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get particle state dimension
    const int statedim = container->GetStateDim(PARTICLEENGINE::Acceleration);

    // get pointer to particle states
    const double* radius = container->GetPtrToState(PARTICLEENGINE::Radius, 0);
    const double* mass = container->GetPtrToState(PARTICLEENGINE::Mass, 0);
    const double* force = container->GetPtrToState(PARTICLEENGINE::Force, 0);
    const double* moment = container->CondGetPtrToState(PARTICLEENGINE::Moment, 0);
    double* acc = container->GetPtrToState(PARTICLEENGINE::Acceleration, 0);
    double* angacc = container->CondGetPtrToState(PARTICLEENGINE::AngularAcceleration, 0);

    // compute acceleration
    for (int i = 0; i < particlestored; ++i)
      UTILS::VecAddScale(&acc[statedim * i], (1.0 / mass[i]), &force[statedim * i]);

    // compute angular acceleration
    if (angacc and moment)
    {
      for (int i = 0; i < particlestored; ++i)
        UTILS::VecAddScale(&angacc[statedim * i],
            (5.0 / (2.0 * mass[i] * UTILS::Pow<2>(radius[i]))), &moment[statedim * i]);
    }
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateParticleEnergy() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateParticleEnergy");

  // evaluate particle kinetic energy contribution
  std::vector<double> kinenergy(1, 0.0);
  {
    std::vector<double> localkinenergy(1, 0.0);
    EvaluateParticleKineticEnergy(localkinenergy[0]);
    comm_.SumAll(localkinenergy.data(), kinenergy.data(), 1);
  }

  // evaluate particle gravitational potential energy contribution
  std::vector<double> gravpotenergy(1, 0.0);
  {
    std::vector<double> localgravpotenergy(1, 0.0);
    EvaluateParticleGravitationalPotentialEnergy(localgravpotenergy[0]);
    comm_.SumAll(localgravpotenergy.data(), gravpotenergy.data(), 1);
  }

  // evaluate elastic potential energy contribution
  std::vector<double> elastpotenergy(1, 0.0);
  {
    std::vector<double> localelastpotenergy(1, 0.0);
    contact_->EvaluateElasticPotentialEnergy(localelastpotenergy[0]);
    comm_.SumAll(localelastpotenergy.data(), elastpotenergy.data(), 1);
  }

  // get specific runtime csv writer
  IO::RuntimeCsvWriter* runtime_csv_writer =
      particleinteractionwriter_->GetSpecificRuntimeCsvWriter("particle-energy");

  // append data vector
  runtime_csv_writer->AppendDataVector("kin_energy", kinenergy);
  runtime_csv_writer->AppendDataVector("grav_pot_energy", gravpotenergy);
  runtime_csv_writer->AppendDataVector("elast_pot_energy", elastpotenergy);
}

void PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateParticleKineticEnergy(
    double& kineticenergy) const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateParticleKineticEnergy");

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get particle state dimension
    const int statedim = container->GetStateDim(PARTICLEENGINE::Position);

    // get pointer to particle states
    const double* radius = container->GetPtrToState(PARTICLEENGINE::Radius, 0);
    const double* mass = container->GetPtrToState(PARTICLEENGINE::Mass, 0);
    const double* vel = container->GetPtrToState(PARTICLEENGINE::Velocity, 0);
    double* angvel = container->CondGetPtrToState(PARTICLEENGINE::AngularVelocity, 0);

    // add translational kinetic energy contribution
    for (int i = 0; i < particlestored; ++i)
      kineticenergy += 0.5 * mass[i] * UTILS::VecDot(&vel[statedim * i], &vel[statedim * i]);

    // add rotational kinetic energy contribution
    if (angvel)
    {
      for (int i = 0; i < particlestored; ++i)
        kineticenergy += 0.5 * (0.4 * mass[i] * UTILS::Pow<2>(radius[i])) *
                         UTILS::VecDot(&angvel[statedim * i], &angvel[statedim * i]);
    }
  }
}

void PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateParticleGravitationalPotentialEnergy(
    double& gravitationalpotentialenergy) const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::ParticleInteractionDEM::EvaluateParticleGravitationalPotentialEnergy");

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get particle state dimension
    const int statedim = container->GetStateDim(PARTICLEENGINE::Position);

    // get pointer to particle states
    const double* pos = container->GetPtrToState(PARTICLEENGINE::Position, 0);
    const double* mass = container->GetPtrToState(PARTICLEENGINE::Mass, 0);

    // add gravitational potential energy contribution
    for (int i = 0; i < particlestored; ++i)
      gravitationalpotentialenergy -= mass[i] * UTILS::VecDot(gravity_.data(), &pos[statedim * i]);
  }
}

BACI_NAMESPACE_CLOSE
