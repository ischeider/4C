/*---------------------------------------------------------------------------*/
/*! \file
\brief smoothed particle hydrodynamics (SPH) interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_interaction_sph.H"

#include "baci_io_pstream.H"
#include "baci_particle_engine_container.H"
#include "baci_particle_engine_interface.H"
#include "baci_particle_interaction_material_handler.H"
#include "baci_particle_interaction_sph_boundary_particle.H"
#include "baci_particle_interaction_sph_density.H"
#include "baci_particle_interaction_sph_equationofstate.H"
#include "baci_particle_interaction_sph_equationofstate_bundle.H"
#include "baci_particle_interaction_sph_kernel.H"
#include "baci_particle_interaction_sph_momentum.H"
#include "baci_particle_interaction_sph_neighbor_pairs.H"
#include "baci_particle_interaction_sph_open_boundary.H"
#include "baci_particle_interaction_sph_phase_change.H"
#include "baci_particle_interaction_sph_pressure.H"
#include "baci_particle_interaction_sph_rigid_particle_contact.H"
#include "baci_particle_interaction_sph_surface_tension.H"
#include "baci_particle_interaction_sph_temperature.H"
#include "baci_particle_interaction_sph_virtual_wall_particle.H"
#include "baci_particle_interaction_utils.H"
#include "baci_particle_wall_interface.H"

#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionSPH::ParticleInteractionSPH(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::ParticleInteractionBase(comm, params), params_sph_(params.sublist("SPH"))
{
  // empty constructor
}

PARTICLEINTERACTION::ParticleInteractionSPH::~ParticleInteractionSPH() = default;

void PARTICLEINTERACTION::ParticleInteractionSPH::Init()
{
  // call base class init
  ParticleInteractionBase::Init();

  // init kernel handler
  InitKernelHandler();

  // init equation of state bundle
  InitEquationOfStateBundle();

  // init neighbor pair handler
  InitNeighborPairHandler();

  // init density handler
  InitDensityHandler();

  // init pressure handler
  InitPressureHandler();

  // init temperature handler
  InitTemperatureHandler();

  // init momentum handler
  InitMomentumHandler();

  // init surface tension handler
  InitSurfaceTensionHandler();

  // init boundary particle handler
  InitBoundaryParticleHandler();

  // init dirichlet open boundary handler
  InitDirichletOpenBoundaryHandler();

  // init neumann open boundary handler
  InitNeumannOpenBoundaryHandler();

  // init virtual wall particle handler
  InitVirtualWallParticleHandler();

  // init phase change handler
  InitPhaseChangeHandler();

  // init rigid particle contact handler
  InitRigidParticleContactHandler();

  // safety check
  if (surfacetension_ and virtualwallparticle_)
    dserror("surface tension formulation with wall interaction not implemented!");
}

void PARTICLEINTERACTION::ParticleInteractionSPH::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // call base class setup
  ParticleInteractionBase::Setup(particleengineinterface, particlewallinterface);

  // setup kernel handler
  kernel_->Setup();

  // setup equation of state bundle
  equationofstatebundle_->Setup();

  // setup neighbor pair handler
  neighborpairs_->Setup(particleengineinterface, particlewallinterface, kernel_);

  // setup density handler
  density_->Setup(particleengineinterface, particlewallinterface, kernel_, particlematerial_,
      equationofstatebundle_, neighborpairs_, virtualwallparticle_);

  // setup pressure handler
  pressure_->Setup(particleengineinterface, particlematerial_, equationofstatebundle_);

  // setup temperature handler
  if (temperature_) temperature_->Setup(particleengineinterface, particlematerial_, neighborpairs_);

  // setup momentum handler
  momentum_->Setup(particleengineinterface, particlewallinterface, kernel_, particlematerial_,
      particleinteractionwriter_, equationofstatebundle_, neighborpairs_, virtualwallparticle_);

  // setup surface tension handler
  if (surfacetension_)
    surfacetension_->Setup(particleengineinterface, kernel_, particlematerial_,
        equationofstatebundle_, neighborpairs_);

  // setup boundary particle handler
  if (boundaryparticle_) boundaryparticle_->Setup(particleengineinterface, neighborpairs_);

  // setup dirichlet open boundary handler
  if (dirichletopenboundary_)
    dirichletopenboundary_->Setup(particleengineinterface, kernel_, particlematerial_,
        equationofstatebundle_, neighborpairs_);

  // setup neumann open boundary handler
  if (neumannopenboundary_)
    neumannopenboundary_->Setup(particleengineinterface, kernel_, particlematerial_,
        equationofstatebundle_, neighborpairs_);

  // setup virtual wall particle handler
  if (virtualwallparticle_)
    virtualwallparticle_->Setup(
        particleengineinterface, particlewallinterface, kernel_, neighborpairs_);

  // setup phase change handler
  if (phasechange_)
    phasechange_->Setup(particleengineinterface, particlematerial_, equationofstatebundle_);

  // setup rigid particle contact handler
  if (rigidparticlecontact_)
    rigidparticlecontact_->Setup(
        particleengineinterface, particlewallinterface, particleinteractionwriter_, neighborpairs_);

  // short screen output
  if ((dirichletopenboundary_ or neumannopenboundary_) and
      particleengineinterface_->HavePeriodicBoundaryConditions())
  {
    if (myrank_ == 0)
      IO::cout << "Warning: periodic boundary and open boundary conditions applied!" << IO::endl;
  }
}

void PARTICLEINTERACTION::ParticleInteractionSPH::WriteRestart() const
{
  // call base class function
  ParticleInteractionBase::WriteRestart();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // call base class function
  ParticleInteractionBase::ReadRestart(reader);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase)
    {
      // insert states of boundary and rigid particles
      particlestates.insert({PARTICLEENGINE::Mass, PARTICLEENGINE::Radius,
          PARTICLEENGINE::BoundaryPressure, PARTICLEENGINE::BoundaryVelocity});
    }
    else if (type == PARTICLEENGINE::DirichletPhase or type == PARTICLEENGINE::NeumannPhase)
    {
      // insert states of open boundary particles
      particlestates.insert({PARTICLEENGINE::Mass, PARTICLEENGINE::Radius, PARTICLEENGINE::Density,
          PARTICLEENGINE::Pressure});
    }
    else
    {
      // insert states of regular phase particles
      particlestates.insert({PARTICLEENGINE::Mass, PARTICLEENGINE::Radius, PARTICLEENGINE::Density,
          PARTICLEENGINE::Pressure});
    }
  }

  // states for density evaluation scheme
  density_->InsertParticleStatesOfParticleTypes(particlestatestotypes);

  // states for temperature evaluation scheme
  if (temperature_) temperature_->InsertParticleStatesOfParticleTypes(particlestatestotypes);

  // insert momentum evaluation dependent states
  momentum_->InsertParticleStatesOfParticleTypes(particlestatestotypes);

  // additional states for surface tension formulation
  if (surfacetension_) surfacetension_->InsertParticleStatesOfParticleTypes(particlestatestotypes);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::SetInitialStates()
{
  // get kernel space dimension
  int kernelspacedim = 0;
  kernel_->KernelSpaceDimension(kernelspacedim);

  // get initial particle spacing
  const double initialparticlespacing = params_sph_.get<double>("INITIALPARTICLESPACING");

  // safety check
  if (not(initialparticlespacing > 0.0)) dserror("negative initial particle spacing!");

  // compute initial particle volume
  const double initialparticlevolume = std::pow(initialparticlespacing, kernelspacedim);

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

    // initial density of current phase
    std::vector<double> initdensity(1);
    initdensity[0] = material->initDensity_;

    // (initial) mass of current phase
    std::vector<double> initmass(1);
    initmass[0] = initdensity[0] * initialparticlevolume;

    // (initial) radius of current phase
    std::vector<double> initradius(1);
    initradius[0] = material->initRadius_;

    // set initial density for respective particles of current type
    if (container->HaveStoredState(PARTICLEENGINE::Density))
      container->SetState(initdensity, PARTICLEENGINE::Density);

    // set initial mass and radius for all particles of current type
    container->SetState(initmass, PARTICLEENGINE::Mass);
    container->SetState(initradius, PARTICLEENGINE::Radius);

    // evaluate initial inertia for respective particles of current type
    if (container->HaveStoredState(PARTICLEENGINE::Inertia))
    {
      // (initial) inertia of current phase
      std::vector<double> initinertia(1);

      if (kernelspacedim == 2)
      {
        // effective particle radius considering initial particle volume in disk shape
        const double effectiveradius = std::sqrt(M_1_PI * initialparticlevolume);

        // inertia for disk shape
        initinertia[0] = 0.5 * initmass[0] * UTILS::Pow<2>(effectiveradius);
      }
      else if (kernelspacedim == 3)
      {
        // effective particle radius considering initial particle volume in spherical shape
        const double effectiveradius = std::pow(0.75 * M_1_PI * initialparticlevolume, 1.0 / 3.0);

        // inertia for spherical shape
        initinertia[0] = 0.4 * initmass[0] * UTILS::Pow<2>(effectiveradius);
      }
      else
      {
        dserror("inertia for particles only in two and three dimensional evaluation given!");
      }

      // set initial inertia for respective particles of current type
      container->SetState(initinertia, PARTICLEENGINE::Inertia);
    }

    // initial states for temperature evaluation
    if (temperature_)
    {
      // get material for current particle type
      const MAT::PAR::ParticleMaterialThermo* material =
          dynamic_cast<const MAT::PAR::ParticleMaterialThermo*>(
              particlematerial_->GetPtrToParticleMatParameter(type_i));

      // initial temperature of current phase
      std::vector<double> inittemperature(1);
      inittemperature[0] = material->initTemperature_;

      // set initial temperature for all particles of current type
      container->SetState(inittemperature, PARTICLEENGINE::Temperature);
    }
  }
}

void PARTICLEINTERACTION::ParticleInteractionSPH::PreEvaluateTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionSPH::PreEvaluateTimeStep");

  // prescribe open boundary states
  if (dirichletopenboundary_) dirichletopenboundary_->PrescribeOpenBoundaryStates(time_);

  // prescribe open boundary states
  if (neumannopenboundary_) neumannopenboundary_->PrescribeOpenBoundaryStates(time_);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::EvaluateInteractions()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionSPH::EvaluateInteractions");

  // evaluate particle neighbor pairs
  neighborpairs_->EvaluateNeighborPairs();

  // init relative positions of virtual particles
  if (virtualwallparticle_)
    virtualwallparticle_->InitRelativePositionsOfVirtualParticles(MaxInteractionDistance());

  // compute density field
  density_->ComputeDensity();

  // compute pressure using equation of state and density
  pressure_->ComputePressure();

  // compute interface quantities
  if (surfacetension_) surfacetension_->ComputeInterfaceQuantities();

  // compute temperature field
  if (temperature_) temperature_->ComputeTemperature();

  // interpolate open boundary states
  if (dirichletopenboundary_) dirichletopenboundary_->InterpolateOpenBoundaryStates();

  // interpolate open boundary states
  if (neumannopenboundary_) neumannopenboundary_->InterpolateOpenBoundaryStates();

  // init boundary particle states
  if (boundaryparticle_) boundaryparticle_->InitBoundaryParticleStates(gravity_);

  // init states at wall contact points
  if (virtualwallparticle_) virtualwallparticle_->InitStatesAtWallContactPoints(gravity_);

  // add momentum contribution to acceleration field
  momentum_->AddAccelerationContribution();

  // add surface tension contribution to acceleration field
  if (surfacetension_) surfacetension_->AddAccelerationContribution();

  // add rigid particle contact contribution to force field
  if (rigidparticlecontact_) rigidparticlecontact_->AddForceContribution();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::PostEvaluateTimeStep(
    std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::ParticleInteractionSPH::PostEvaluateTimeStep");

  // check open boundary phase change
  if (dirichletopenboundary_)
    dirichletopenboundary_->CheckOpenBoundaryPhaseChange(MaxInteractionDistance());

  // check open boundary phase change
  if (neumannopenboundary_)
    neumannopenboundary_->CheckOpenBoundaryPhaseChange(MaxInteractionDistance());

  // evaluate phase change
  if (phasechange_) phasechange_->EvaluatePhaseChange(particlesfromphasetophase);
}

double PARTICLEINTERACTION::ParticleInteractionSPH::MaxInteractionDistance() const
{
  return MaxParticleRadius();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::DistributeInteractionHistory() const
{
  // nothing to do
}

void PARTICLEINTERACTION::ParticleInteractionSPH::CommunicateInteractionHistory() const
{
  // nothing to do
}

void PARTICLEINTERACTION::ParticleInteractionSPH::SetCurrentTime(const double currenttime)
{
  // call base class method
  ParticleInteractionBase::SetCurrentTime(currenttime);

  // set current time
  if (temperature_) temperature_->SetCurrentTime(currenttime);

  // set current time
  if (surfacetension_) surfacetension_->SetCurrentTime(currenttime);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::SetCurrentStepSize(const double currentstepsize)
{
  // call base class method
  ParticleInteractionBase::SetCurrentStepSize(currentstepsize);

  // set current step size
  density_->SetCurrentStepSize(currentstepsize);

  // set current step size
  if (temperature_) temperature_->SetCurrentStepSize(currentstepsize);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitKernelHandler()
{
  // get type of smoothed particle hydrodynamics kernel
  INPAR::PARTICLE::KernelType kerneltype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::KernelType>(params_sph_, "KERNEL");

  // create kernel handler
  switch (kerneltype)
  {
    case INPAR::PARTICLE::CubicSpline:
    {
      kernel_ = std::make_shared<PARTICLEINTERACTION::SPHKernelCubicSpline>(params_sph_);
      break;
    }
    case INPAR::PARTICLE::QuinticSpline:
    {
      kernel_ = std::make_shared<PARTICLEINTERACTION::SPHKernelQuinticSpline>(params_sph_);
      break;
    }
    default:
    {
      dserror("unknown kernel type!");
      break;
    }
  }

  // init kernel handler
  kernel_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitEquationOfStateBundle()
{
  // create equation of state bundle
  equationofstatebundle_ =
      std::make_shared<PARTICLEINTERACTION::SPHEquationOfStateBundle>(params_sph_);

  // init equation of state bundle
  equationofstatebundle_->Init(particlematerial_);
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitNeighborPairHandler()
{
  // create neighbor pair handler
  neighborpairs_ = std::make_shared<PARTICLEINTERACTION::SPHNeighborPairs>();

  // init neighbor pair handler
  neighborpairs_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitDensityHandler()
{
  // get type of smoothed particle hydrodynamics density evaluation scheme
  INPAR::PARTICLE::DensityEvaluationScheme densityevaluationscheme =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DensityEvaluationScheme>(
          params_sph_, "DENSITYEVALUATION");

  // create density handler
  switch (densityevaluationscheme)
  {
    case INPAR::PARTICLE::DensitySummation:
    {
      density_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensitySummation>(
          new PARTICLEINTERACTION::SPHDensitySummation(params_sph_));
      break;
    }
    case INPAR::PARTICLE::DensityIntegration:
    {
      density_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityIntegration>(
          new PARTICLEINTERACTION::SPHDensityIntegration(params_sph_));
      break;
    }
    case INPAR::PARTICLE::DensityPredictCorrect:
    {
      density_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityPredictCorrect>(
          new PARTICLEINTERACTION::SPHDensityPredictCorrect(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown density evaluation scheme type!");
      break;
    }
  }

  // init density handler
  density_->Init();

  // safety check
  if (densityevaluationscheme != INPAR::PARTICLE::DensityPredictCorrect and
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DensityCorrectionScheme>(
          params_sph_, "DENSITYCORRECTION") != INPAR::PARTICLE::NoCorrection)
    dserror(
        "the density correction scheme set is not valid with the current density evaluation "
        "scheme!");
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitPressureHandler()
{
  // create pressure handler
  pressure_ =
      std::unique_ptr<PARTICLEINTERACTION::SPHPressure>(new PARTICLEINTERACTION::SPHPressure());

  // init pressure handler
  pressure_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitTemperatureHandler()
{
  // get type of smoothed particle hydrodynamics temperature evaluation scheme
  INPAR::PARTICLE::TemperatureEvaluationScheme temperatureevaluationscheme =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::TemperatureEvaluationScheme>(
          params_sph_, "TEMPERATUREEVALUATION");

  // create temperature handler
  switch (temperatureevaluationscheme)
  {
    case INPAR::PARTICLE::NoTemperatureEvaluation:
    {
      temperature_ = std::unique_ptr<PARTICLEINTERACTION::SPHTemperature>(nullptr);
      break;
    }
    case INPAR::PARTICLE::TemperatureIntegration:
    {
      temperature_ = std::unique_ptr<PARTICLEINTERACTION::SPHTemperature>(
          new PARTICLEINTERACTION::SPHTemperature(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown temperature evaluation scheme!");
      break;
    }
  }

  // init temperature handler
  if (temperature_) temperature_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitMomentumHandler()
{
  // create momentum handler
  momentum_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentum>(
      new PARTICLEINTERACTION::SPHMomentum(params_sph_));

  // init momentum handler
  momentum_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitSurfaceTensionHandler()
{
  // get type of smoothed particle hydrodynamics surface tension formulation
  INPAR::PARTICLE::SurfaceTensionFormulation surfacetensionformulation =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::SurfaceTensionFormulation>(
          params_sph_, "SURFACETENSIONFORMULATION");

  // create surface tension handler
  switch (surfacetensionformulation)
  {
    case INPAR::PARTICLE::NoSurfaceTension:
    {
      surfacetension_ = std::unique_ptr<PARTICLEINTERACTION::SPHSurfaceTension>(nullptr);
      break;
    }
    case INPAR::PARTICLE::ContinuumSurfaceForce:
    {
      surfacetension_ = std::unique_ptr<PARTICLEINTERACTION::SPHSurfaceTension>(
          new PARTICLEINTERACTION::SPHSurfaceTension(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown surface tension formulation type!");
      break;
    }
  }

  // init surface tension handler
  if (surfacetension_) surfacetension_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitBoundaryParticleHandler()
{
  // get type of boundary particle formulation
  INPAR::PARTICLE::BoundaryParticleFormulationType boundaryparticleformulation =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::BoundaryParticleFormulationType>(
          params_sph_, "BOUNDARYPARTICLEFORMULATION");

  // create boundary particle handler
  switch (boundaryparticleformulation)
  {
    case INPAR::PARTICLE::NoBoundaryFormulation:
    {
      boundaryparticle_ = std::unique_ptr<PARTICLEINTERACTION::SPHBoundaryParticleBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::AdamiBoundaryFormulation:
    {
      boundaryparticle_ = std::unique_ptr<PARTICLEINTERACTION::SPHBoundaryParticleAdami>(
          new PARTICLEINTERACTION::SPHBoundaryParticleAdami(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown boundary particle formulation type!");
      break;
    }
  }

  // init boundary particle handler
  if (boundaryparticle_) boundaryparticle_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitDirichletOpenBoundaryHandler()
{
  // get type of dirichlet open boundary
  INPAR::PARTICLE::DirichletOpenBoundaryType dirichletopenboundarytype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DirichletOpenBoundaryType>(
          params_sph_, "DIRICHLETBOUNDARYTYPE");

  // create open boundary handler
  switch (dirichletopenboundarytype)
  {
    case INPAR::PARTICLE::NoDirichletOpenBoundary:
    {
      dirichletopenboundary_ = std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::DirichletNormalToPlane:
    {
      dirichletopenboundary_ = std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryDirichlet>(
          new PARTICLEINTERACTION::SPHOpenBoundaryDirichlet(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown dirichlet open boundary type!");
      break;
    }
  }

  // init open boundary handler
  if (dirichletopenboundary_) dirichletopenboundary_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitNeumannOpenBoundaryHandler()
{
  // get type of neumann open boundary
  INPAR::PARTICLE::NeumannOpenBoundaryType neumannopenboundarytype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::NeumannOpenBoundaryType>(
          params_sph_, "NEUMANNBOUNDARYTYPE");

  // create open boundary handler
  switch (neumannopenboundarytype)
  {
    case INPAR::PARTICLE::NoNeumannOpenBoundary:
    {
      neumannopenboundary_ = std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::NeumannNormalToPlane:
    {
      neumannopenboundary_ = std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryNeumann>(
          new PARTICLEINTERACTION::SPHOpenBoundaryNeumann(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown neumann open boundary type!");
      break;
    }
  }

  // init open boundary handler
  if (neumannopenboundary_) neumannopenboundary_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitVirtualWallParticleHandler()
{
  // get type of wall formulation
  INPAR::PARTICLE::WallFormulationType wallformulation =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallFormulationType>(
          params_sph_, "WALLFORMULATION");

  // create virtual wall particle handler
  switch (wallformulation)
  {
    case INPAR::PARTICLE::NoWallFormulation:
    {
      virtualwallparticle_ = std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle>(nullptr);
      break;
    }
    case INPAR::PARTICLE::VirtualParticleWallFormulation:
    {
      virtualwallparticle_ =
          std::make_shared<PARTICLEINTERACTION::SPHVirtualWallParticle>(params_sph_);
      break;
    }
    default:
    {
      dserror("unknown wall formulation type!");
      break;
    }
  }

  // init virtual wall particle handler
  if (virtualwallparticle_) virtualwallparticle_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitPhaseChangeHandler()
{
  // get type of phase change
  INPAR::PARTICLE::PhaseChangeType phasechangetype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::PhaseChangeType>(params_sph_, "PHASECHANGETYPE");

  // create phase change handler
  switch (phasechangetype)
  {
    case INPAR::PARTICLE::NoPhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::OneWayScalarBelowToAbovePhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarBelowToAbove>(
          new PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarBelowToAbove(params_sph_));
      break;
    }
    case INPAR::PARTICLE::OneWayScalarAboveToBelowPhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarAboveToBelow>(
          new PARTICLEINTERACTION::SPHPhaseChangeOneWayScalarAboveToBelow(params_sph_));
      break;
    }
    case INPAR::PARTICLE::TwoWayScalarPhaseChange:
    {
      phasechange_ = std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar>(
          new PARTICLEINTERACTION::SPHPhaseChangeTwoWayScalar(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown boundary particle formulation type!");
      break;
    }
  }

  // init phase change handler
  if (phasechange_) phasechange_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionSPH::InitRigidParticleContactHandler()
{
  // get type of rigid particle contact
  INPAR::PARTICLE::RigidParticleContactType rigidparticlecontacttype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::RigidParticleContactType>(
          params_sph_, "RIGIDPARTICLECONTACTTYPE");

  // create rigid particle contact handler
  switch (rigidparticlecontacttype)
  {
    case INPAR::PARTICLE::NoRigidParticleContact:
    {
      rigidparticlecontact_ =
          std::unique_ptr<PARTICLEINTERACTION::SPHRigidParticleContactBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::ElasticRigidParticleContact:
    {
      rigidparticlecontact_ = std::unique_ptr<PARTICLEINTERACTION::SPHRigidParticleContactElastic>(
          new PARTICLEINTERACTION::SPHRigidParticleContactElastic(params_sph_));
      break;
    }
    default:
    {
      dserror("unknown rigid particle contact type!");
      break;
    }
  }

  // init rigid particle contact handler
  if (rigidparticlecontact_) rigidparticlecontact_->Init();
}

BACI_NAMESPACE_CLOSE
