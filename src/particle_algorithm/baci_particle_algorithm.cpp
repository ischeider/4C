/*---------------------------------------------------------------------------*/
/*! \file
\brief algorithm to control particle simulations
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_algorithm.H"

#include "baci_inpar_particle.H"
#include "baci_io.H"
#include "baci_io_pstream.H"
#include "baci_lib_resulttest.H"
#include "baci_particle_algorithm_gravity.H"
#include "baci_particle_algorithm_growth.H"
#include "baci_particle_algorithm_initial_field.H"
#include "baci_particle_algorithm_input_generator.H"
#include "baci_particle_algorithm_result_test.H"
#include "baci_particle_algorithm_timint.H"
#include "baci_particle_algorithm_utils.H"
#include "baci_particle_algorithm_viscous_damping.H"
#include "baci_particle_engine.H"
#include "baci_particle_engine_communication_utils.H"
#include "baci_particle_engine_container.H"
#include "baci_particle_engine_object.H"
#include "baci_particle_interaction_base.H"
#include "baci_particle_interaction_dem.H"
#include "baci_particle_interaction_sph.H"
#include "baci_particle_rigidbody.H"
#include "baci_particle_rigidbody_result_test.H"
#include "baci_particle_wall.H"
#include "baci_particle_wall_result_test.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::ParticleAlgorithm::ParticleAlgorithm(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : AlgorithmBase(comm, params),
      myrank_(comm.MyPID()),
      params_(params),
      numparticlesafterlastloadbalance_(0),
      transferevery_(DRT::INPUT::IntegralValue<int>(params_, "TRANSFER_EVERY")),
      writeresultsevery_(params.get<int>("RESULTSEVRY")),
      writerestartevery_(params.get<int>("RESTARTEVRY")),
      writeresultsthisstep_(true),
      writerestartthisstep_(false),
      isrestarted_(false)
{
  // empty constructor
}

PARTICLEALGORITHM::ParticleAlgorithm::~ParticleAlgorithm() = default;

void PARTICLEALGORITHM::ParticleAlgorithm::Init(
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles)
{
  // init particle engine
  InitParticleEngine();

  // init particle wall handler
  InitParticleWall();

  // init rigid body handler
  InitParticleRigidBody();

  // init particle time integration
  InitParticleTimeIntegration();

  // init particle interaction handler
  InitParticleInteraction();

  // init particle gravity handler
  InitParticleGravity();

  // init particle growth handler
  InitParticleGrowth();

  // init viscous damping handler
  InitViscousDamping();

  // set initial particles to vector of particles to be distributed
  particlestodistribute_ = initialparticles;

  // clear vector of initial particles in global problem
  initialparticles.clear();
}

void PARTICLEALGORITHM::ParticleAlgorithm::Setup()
{
  // generate initial particles
  if (not isrestarted_) GenerateInitialParticles();

  // determine all particle types
  DetermineParticleTypes();

  // determine particle states of all particle types
  DetermineParticleStatesOfParticleTypes();

  // setup particle engine
  particleengine_->Setup(particlestatestotypes_);

  // setup wall handler
  if (particlewall_) particlewall_->Setup(particleengine_);

  // setup rigid body handler
  if (particlerigidbody_) particlerigidbody_->Setup(particleengine_);

  // setup particle time integration
  particletimint_->Setup(particleengine_, particlerigidbody_);

  // setup particle interaction handler
  if (particleinteraction_) particleinteraction_->Setup(particleengine_, particlewall_);

  // setup gravity handler
  if (particlegravity_) particlegravity_->Setup();

  // setup particle growth handler
  if (particlegrowth_) particlegrowth_->Setup(particleengine_);

  // setup viscous damping handler
  if (viscousdamping_) viscousdamping_->Setup(particleengine_);

  // setup initial particles
  SetupInitialParticles();

  // setup initial rigid bodies
  if (particlerigidbody_) SetupInitialRigidBodies();

  // distribute load among processors
  DistributeLoadAmongProcs();

  // ghost particles on other processors
  particleengine_->GhostParticles();

  // build global id to local index map
  particleengine_->BuildGlobalIDToLocalIndexMap();

  // build potential neighbor relation
  if (particleinteraction_) BuildPotentialNeighborRelation();

  // setup initial states
  if (not isrestarted_) SetupInitialStates();

  // write initial output
  if (not isrestarted_) WriteOutput();
}

void PARTICLEALGORITHM::ParticleAlgorithm::ReadRestart(const int restartstep)
{
  // clear vector of particles to be distributed
  particlestodistribute_.clear();

  // create discretization reader
  const std::shared_ptr<IO::DiscretizationReader> reader =
      particleengine_->BinDisReader(restartstep);

  // safety check
  if (restartstep != reader->ReadInt("step")) dserror("time step on file not equal to given step!");

  // get restart time
  double restarttime = reader->ReadDouble("time");

  // read restart of particle engine
  particleengine_->ReadRestart(reader, particlestodistribute_);

  // read restart of rigid body handler
  if (particlerigidbody_) particlerigidbody_->ReadRestart(reader);

  // read restart of particle interaction handler
  if (particleinteraction_) particleinteraction_->ReadRestart(reader);

  // read restart of wall handler
  if (particlewall_) particlewall_->ReadRestart(restartstep);

  // set time and step after restart
  SetTimeStep(restarttime, restartstep);

  // set flag indicating restart to true
  isrestarted_ = true;

  // short screen output
  if (myrank_ == 0)
    IO::cout << "====== restart of the particle simulation from step " << restartstep << IO::endl;
}

void PARTICLEALGORITHM::ParticleAlgorithm::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    PrepareTimeStep();

    // pre evaluate time step
    PreEvaluateTimeStep();

    // integrate time step
    IntegrateTimeStep();

    // post evaluate time step
    PostEvaluateTimeStep();

    // write output
    WriteOutput();

    // write restart information
    WriteRestart();
  }
}

void PARTICLEALGORITHM::ParticleAlgorithm::PrepareTimeStep(bool print_header)
{
  // increment time and step
  IncrementTimeAndStep();

  // set current time
  SetCurrentTime();

  // set current step size
  SetCurrentStepSize();

  // print header
  if (print_header) PrintHeader();

  // update result and restart control flags
  writeresultsthisstep_ = (writeresultsevery_ and (Step() % writeresultsevery_ == 0));
  writerestartthisstep_ = (writerestartevery_ and (Step() % writerestartevery_ == 0));

  // set current write result flag
  SetCurrentWriteResultFlag();
}

void PARTICLEALGORITHM::ParticleAlgorithm::PreEvaluateTimeStep()
{
  // pre evaluate time step
  if (particleinteraction_) particleinteraction_->PreEvaluateTimeStep();

  if (particlegrowth_) particlegrowth_->PreEvaluateTimeStep(Time());
}

void PARTICLEALGORITHM::ParticleAlgorithm::IntegrateTimeStep()
{
  // time integration scheme specific pre-interaction routine
  particletimint_->PreInteractionRoutine();

  // update connectivity
  UpdateConnectivity();

  // evaluate time step
  EvaluateTimeStep();

  // time integration scheme specific post-interaction routine
  particletimint_->PostInteractionRoutine();
}

void PARTICLEALGORITHM::ParticleAlgorithm::PostEvaluateTimeStep()
{
  // post evaluate time step
  std::vector<PARTICLEENGINE::ParticleTypeToType> particlesfromphasetophase;
  if (particleinteraction_) particleinteraction_->PostEvaluateTimeStep(particlesfromphasetophase);

  if (particlerigidbody_)
  {
    // have rigid body phase change
    if (particlerigidbody_->HaveRigidBodyPhaseChange(particlesfromphasetophase))
    {
      // update connectivity
      UpdateConnectivity();

      // evaluate rigid body phase change
      particlerigidbody_->EvaluateRigidBodyPhaseChange(particlesfromphasetophase);
    }
  }
}

void PARTICLEALGORITHM::ParticleAlgorithm::WriteOutput() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::ParticleAlgorithm::WriteOutput");

  // write result step
  if (writeresultsthisstep_)
  {
    // write particle runtime output
    particleengine_->WriteParticleRuntimeOutput(Step(), Time());

    // write binning discretization output (debug feature)
    particleengine_->WriteBinDisOutput(Step(), Time());

    // write rigid body runtime output
    if (particlerigidbody_) particlerigidbody_->WriteRigidBodyRuntimeOutput(Step(), Time());

    // write interaction runtime output
    if (particleinteraction_) particleinteraction_->WriteInteractionRuntimeOutput(Step(), Time());

    // write wall runtime output
    if (particlewall_) particlewall_->WriteWallRuntimeOutput(Step(), Time());
  }
}

void PARTICLEALGORITHM::ParticleAlgorithm::WriteRestart() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::ParticleAlgorithm::WriteRestart");

  // write restart step
  if (writerestartthisstep_)
  {
    // write restart of particle engine
    particleengine_->WriteRestart(Step(), Time());

    // write restart of rigid body handler
    if (particlerigidbody_) particlerigidbody_->WriteRestart();

    // write restart of particle interaction handler
    if (particleinteraction_) particleinteraction_->WriteRestart();

    // write restart of wall handler
    if (particlewall_) particlewall_->WriteRestart(Step(), Time());

    // short screen output
    if (myrank_ == 0)
      IO::cout(IO::verbose) << "====== restart of the particle simulation written in step "
                            << Step() << IO::endl;
  }
}

std::vector<std::shared_ptr<DRT::ResultTest>>
PARTICLEALGORITHM::ParticleAlgorithm::CreateResultTests()
{
  // build global id to local index map
  particleengine_->BuildGlobalIDToLocalIndexMap();

  // particle field specific result test objects
  std::vector<std::shared_ptr<DRT::ResultTest>> allresulttests(0);

  // particle result test
  {
    // create and init particle result test
    std::shared_ptr<PARTICLEALGORITHM::ParticleResultTest> particleresulttest =
        std::make_shared<PARTICLEALGORITHM::ParticleResultTest>();
    particleresulttest->Init();

    // setup particle result test
    particleresulttest->Setup(particleengine_);

    allresulttests.push_back(particleresulttest);
  }

  // wall result test
  if (particlewall_)
  {
    // create and init wall result test
    std::shared_ptr<PARTICLEWALL::WallResultTest> wallresulttest =
        std::make_shared<PARTICLEWALL::WallResultTest>();
    wallresulttest->Init();

    // setup wall result test
    wallresulttest->Setup(particlewall_);

    allresulttests.push_back(wallresulttest);
  }

  if (particlerigidbody_)
  {
    // create and init rigid body result test
    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyResultTest> rigidbodyresulttest =
        std::make_shared<PARTICLERIGIDBODY::RigidBodyResultTest>();
    rigidbodyresulttest->Init();

    // setup rigid body result test
    rigidbodyresulttest->Setup(particlerigidbody_);

    allresulttests.push_back(rigidbodyresulttest);
  }

  return allresulttests;
}

void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleEngine()
{
  // create and init particle engine
  particleengine_ = std::make_shared<PARTICLEENGINE::ParticleEngine>(Comm(), params_);
  particleengine_->Init();
}

void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleWall()
{
  // get type of particle wall source
  INPAR::PARTICLE::ParticleWallSource particlewallsource =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleWallSource>(
          params_, "PARTICLE_WALL_SOURCE");

  // create particle wall handler
  switch (particlewallsource)
  {
    case INPAR::PARTICLE::NoParticleWall:
    {
      particlewall_ = std::shared_ptr<PARTICLEWALL::WallHandlerBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::DiscretCondition:
    {
      particlewall_ = std::make_shared<PARTICLEWALL::WallHandlerDiscretCondition>(Comm(), params_);
      break;
    }
    case INPAR::PARTICLE::BoundingBox:
    {
      particlewall_ = std::make_shared<PARTICLEWALL::WallHandlerBoundingBox>(Comm(), params_);
      break;
    }
    default:
    {
      dserror("unknown type of particle wall source!");
      break;
    }
  }

  // init particle wall handler
  if (particlewall_) particlewall_->Init(particleengine_->GetBinningStrategy());
}

void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleRigidBody()
{
  // create rigid body handler
  if (DRT::INPUT::IntegralValue<int>(params_, "RIGID_BODY_MOTION"))
    particlerigidbody_ = std::make_shared<PARTICLERIGIDBODY::RigidBodyHandler>(Comm(), params_);

  // init rigid body handler
  if (particlerigidbody_) particlerigidbody_->Init();
}

void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleTimeIntegration()
{
  // get particle time integration scheme
  INPAR::PARTICLE::DynamicType timinttype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DynamicType>(params_, "DYNAMICTYP");

  // create particle time integration
  switch (timinttype)
  {
    case INPAR::PARTICLE::dyna_semiimpliciteuler:
    {
      particletimint_ = std::unique_ptr<PARTICLEALGORITHM::TimIntSemiImplicitEuler>(
          new PARTICLEALGORITHM::TimIntSemiImplicitEuler(params_));
      break;
    }
    case INPAR::PARTICLE::dyna_velocityverlet:
    {
      particletimint_ = std::unique_ptr<PARTICLEALGORITHM::TimIntVelocityVerlet>(
          new PARTICLEALGORITHM::TimIntVelocityVerlet(params_));
      break;
    }
    default:
    {
      dserror("unknown particle time integration scheme!");
      break;
    }
  }

  // init particle time integration
  particletimint_->Init();
}

void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleInteraction()
{
  // get particle interaction type
  INPAR::PARTICLE::InteractionType interactiontype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::InteractionType>(params_, "INTERACTION");

  // create particle interaction handler
  switch (interactiontype)
  {
    case INPAR::PARTICLE::interaction_none:
    {
      particleinteraction_ = std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::interaction_sph:
    {
      particleinteraction_ = std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionSPH>(
          new PARTICLEINTERACTION::ParticleInteractionSPH(Comm(), params_));
      break;
    }
    case INPAR::PARTICLE::interaction_dem:
    {
      particleinteraction_ = std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionDEM>(
          new PARTICLEINTERACTION::ParticleInteractionDEM(Comm(), params_));
      break;
    }
    default:
    {
      dserror("unknown particle interaction type!");
      break;
    }
  }

  // init particle interaction handler
  if (particleinteraction_) particleinteraction_->Init();
}

void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleGravity()
{
  // init gravity acceleration vector
  std::vector<double> gravity;
  std::string value;
  std::istringstream gravitystream(
      Teuchos::getNumericStringParameter(params_, "GRAVITY_ACCELERATION"));

  while (gravitystream >> value) gravity.push_back(std::atof(value.c_str()));

  // safety check
  if (static_cast<int>(gravity.size()) != 3)
    dserror("dimension (dim = %d) of gravity acceleration vector is wrong!",
        static_cast<int>(gravity.size()));

  // get magnitude of gravity
  double temp = 0.0;
  for (double g : gravity) temp += g * g;
  const double gravity_norm = std::sqrt(temp);

  // create particle gravity handler
  if (gravity_norm > 0.0)
    particlegravity_ = std::unique_ptr<PARTICLEALGORITHM::GravityHandler>(
        new PARTICLEALGORITHM::GravityHandler(params_));

  // init particle gravity handler
  if (particlegravity_) particlegravity_->Init(gravity);
}

void PARTICLEALGORITHM::ParticleAlgorithm::InitParticleGrowth()
{
  // get particle growth function
  const int growthfunc = params_.get<int>("RADIUS_RAMP_FUNCT");

  // create particle growth handler
  if (growthfunc > 0)
    particlegrowth_ = std::unique_ptr<PARTICLEALGORITHM::ParticleGrowthHandler>(
        new PARTICLEALGORITHM::ParticleGrowthHandler(growthfunc));

  // init particle growth handler
  if (particlegrowth_) particlegrowth_->Init();
}

void PARTICLEALGORITHM::ParticleAlgorithm::InitViscousDamping()
{
  // get viscous damping factor
  const double viscdampfac = params_.get<double>("VISCOUS_DAMPING");

  // create viscous damping handler
  if (viscdampfac > 0.0)
    viscousdamping_ = std::unique_ptr<PARTICLEALGORITHM::ViscousDampingHandler>(
        new PARTICLEALGORITHM::ViscousDampingHandler(viscdampfac));

  // init viscous damping handler
  if (viscousdamping_) viscousdamping_->Init();
}

void PARTICLEALGORITHM::ParticleAlgorithm::GenerateInitialParticles()
{
  // create and init particle input generator
  std::unique_ptr<PARTICLEALGORITHM::InputGenerator> particleinputgenerator =
      std::unique_ptr<PARTICLEALGORITHM::InputGenerator>(
          new PARTICLEALGORITHM::InputGenerator(Comm(), params_));
  particleinputgenerator->Init();

  // generate particles
  particleinputgenerator->GenerateParticles(particlestodistribute_);
}

void PARTICLEALGORITHM::ParticleAlgorithm::DetermineParticleTypes()
{
  // init map relating particle types to dynamic load balance factor
  std::map<PARTICLEENGINE::TypeEnum, double> typetodynloadbal;

  // read parameters relating particle types to values
  PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues(
      params_, "PHASE_TO_DYNLOADBALFAC", typetodynloadbal);

  // insert into map of particle types and corresponding states with empty set
  for (auto& typeIt : typetodynloadbal)
    particlestatestotypes_.insert(
        std::make_pair(typeIt.first, std::set<PARTICLEENGINE::StateEnum>()));

  // safety check
  for (auto& particle : particlestodistribute_)
    if (not particlestatestotypes_.count(particle->ReturnParticleType()))
      dserror("particle type '%s' of initial particle not defined!",
          PARTICLEENGINE::EnumToTypeName(particle->ReturnParticleType()).c_str());
}

void PARTICLEALGORITHM::ParticleAlgorithm::DetermineParticleStatesOfParticleTypes()
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes_)
  {
    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // insert default particle states
    particlestates.insert({PARTICLEENGINE::Position, PARTICLEENGINE::Velocity,
        PARTICLEENGINE::Acceleration, PARTICLEENGINE::LastTransferPosition});
  }

  // insert integration dependent states of all particle types
  particletimint_->InsertParticleStatesOfParticleTypes(particlestatestotypes_);

  // insert interaction dependent states of all particle types
  if (particleinteraction_)
    particleinteraction_->InsertParticleStatesOfParticleTypes(particlestatestotypes_);

  // insert wall handler dependent states of all particle types
  if (particlewall_) particlewall_->InsertParticleStatesOfParticleTypes(particlestatestotypes_);

  // insert rigid body handler dependent states of all particle types
  if (particlerigidbody_)
    particlerigidbody_->InsertParticleStatesOfParticleTypes(particlestatestotypes_);

  // insert particlegrowth handler dependent state of all particles with growth
  if (particlegrowth_) particlegrowth_->InsertParticleStatesOfParticleTypes(particlestatestotypes_);
}

void PARTICLEALGORITHM::ParticleAlgorithm::SetupInitialParticles()
{
  // get unique global ids for all particles
  if (not isrestarted_) particleengine_->GetUniqueGlobalIdsForAllParticles(particlestodistribute_);

  // erase particles outside bounding box
  particleengine_->EraseParticlesOutsideBoundingBox(particlestodistribute_);

  // distribute particles to owning processor
  particleengine_->DistributeParticles(particlestodistribute_);

  // distribute interaction history
  if (particleinteraction_) particleinteraction_->DistributeInteractionHistory();
}

void PARTICLEALGORITHM::ParticleAlgorithm::SetupInitialRigidBodies()
{
  // set initial affiliation pair data
  if (not isrestarted_) particlerigidbody_->SetInitialAffiliationPairData();

  // set unique global ids for all rigid bodies
  if (not isrestarted_) particlerigidbody_->SetUniqueGlobalIdsForAllRigidBodies();

  // allocate rigid body states
  if (not isrestarted_) particlerigidbody_->AllocateRigidBodyStates();

  // distribute rigid body
  particlerigidbody_->DistributeRigidBody();
}

void PARTICLEALGORITHM::ParticleAlgorithm::SetupInitialStates()
{
  // set initial states
  if (particleinteraction_) particleinteraction_->SetInitialStates();

  // initialize radius information if particle growth
  if (particlegrowth_) particlegrowth_->SetInitialStates();

  // initialize rigid body mass quantities and orientation
  if (particlerigidbody_) particlerigidbody_->InitializeRigidBodyMassQuantitiesAndOrientation();

  // set initial conditions
  SetInitialConditions();

  // time integration scheme specific initialization routine
  particletimint_->SetInitialStates();

  // evaluate consistent initial states
  {
    // pre evaluate time step
    PreEvaluateTimeStep();

    // update connectivity
    UpdateConnectivity();

    // evaluate time step
    EvaluateTimeStep();

    // post evaluate time step
    PostEvaluateTimeStep();
  }
}

void PARTICLEALGORITHM::ParticleAlgorithm::UpdateConnectivity()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::ParticleAlgorithm::UpdateConnectivity");

#ifdef DEBUG
  // check number of unique global ids
  particleengine_->CheckNumberOfUniqueGlobalIds();
#endif

  // check particle interaction distance concerning bin size
  if (particleinteraction_)
    particleinteraction_->CheckParticleInteractionDistanceConcerningBinSize();

  // check that wall nodes are located in bounding box
  if (particlewall_) particlewall_->CheckWallNodesLocatedInBoundingBox();

  if (CheckLoadTransferNeeded())
  {
    // transfer load between processors
    TransferLoadBetweenProcs();

    // distribute load among processors
    if (CheckLoadRedistributionNeeded()) DistributeLoadAmongProcs();

    // ghost particles on other processors
    particleengine_->GhostParticles();

    // build global id to local index map
    particleengine_->BuildGlobalIDToLocalIndexMap();

    // build potential neighbor relation
    if (particleinteraction_) BuildPotentialNeighborRelation();
  }
  else
  {
    // refresh particles being ghosted on other processors
    particleengine_->RefreshParticles();
  }
}

bool PARTICLEALGORITHM::ParticleAlgorithm::CheckLoadTransferNeeded()
{
  bool transferload = transferevery_ or writeresultsthisstep_ or writerestartthisstep_;

  // check max position increment
  transferload |= CheckMaxPositionIncrement();

  // check for valid particle connectivity
  transferload |= (not particleengine_->HaveValidParticleConnectivity());

  // check for valid particle neighbors
  if (particleinteraction_) transferload |= (not particleengine_->HaveValidParticleNeighbors());

  // check for valid wall neighbors
  if (particleinteraction_ and particlewall_)
    transferload |= (not particlewall_->HaveValidWallNeighbors());

  return transferload;
}

bool PARTICLEALGORITHM::ParticleAlgorithm::CheckMaxPositionIncrement()
{
  // get maximum particle interaction distance
  double allprocmaxinteractiondistance = 0.0;
  if (particleinteraction_)
  {
    double maxinteractiondistance = particleinteraction_->MaxInteractionDistance();
    Comm().MaxAll(&maxinteractiondistance, &allprocmaxinteractiondistance, 1);
  }

  // get max particle position increment since last transfer
  double maxparticlepositionincrement = 0.0;
  GetMaxParticlePositionIncrement(maxparticlepositionincrement);

  // get max wall position increment since last transfer
  double maxwallpositionincrement = 0.0;
  if (particlewall_) particlewall_->GetMaxWallPositionIncrement(maxwallpositionincrement);

  // get max overall position increment since last transfer
  const double maxpositionincrement =
      std::max(maxparticlepositionincrement, maxwallpositionincrement);

  // get allowed position increment
  const double allowedpositionincrement =
      0.5 * (particleengine_->MinBinSize() - allprocmaxinteractiondistance);

  // check if a particle transfer is needed based on a worst case scenario:
  // two particles approach each other with maximum position increment in one spatial dimension
  return (maxpositionincrement > allowedpositionincrement);
}

void PARTICLEALGORITHM::ParticleAlgorithm::GetMaxParticlePositionIncrement(
    double& allprocmaxpositionincrement)
{
  // maximum position increment since last particle transfer
  double maxpositionincrement = 0.0;

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengine_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored == 0) continue;

    // get particle state dimension
    int statedim = container->GetStateDim(PARTICLEENGINE::Position);

    // position increment of particle
    double positionincrement[3];

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      // get pointer to particle states
      const double* pos = container->GetPtrToState(PARTICLEENGINE::Position, i);
      const double* lasttransferpos =
          container->GetPtrToState(PARTICLEENGINE::LastTransferPosition, i);

      // position increment of particle considering periodic boundaries
      particleengine_->DistanceBetweenParticles(pos, lasttransferpos, positionincrement);

      // iterate over spatial dimension
      for (int dim = 0; dim < statedim; ++dim)
        maxpositionincrement = std::max(maxpositionincrement, std::abs(positionincrement[dim]));
    }
  }

  // bin size safety check
  if (maxpositionincrement > particleengine_->MinBinSize())
    dserror("a particle traveled more than one bin on this processor!");

  // get maximum particle position increment on all processors
  Comm().MaxAll(&maxpositionincrement, &allprocmaxpositionincrement, 1);
}

void PARTICLEALGORITHM::ParticleAlgorithm::TransferLoadBetweenProcs()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::ParticleAlgorithm::TransferLoadBetweenProcs");

  // transfer particles to new bins and processors
  particleengine_->TransferParticles();

  // transfer wall elements and nodes
  if (particlewall_) particlewall_->TransferWallElementsAndNodes();

  // communicate rigid body
  if (particlerigidbody_) particlerigidbody_->CommunicateRigidBody();

  // communicate interaction history
  if (particleinteraction_) particleinteraction_->CommunicateInteractionHistory();

  // short screen output
  if (myrank_ == 0) IO::cout(IO::verbose) << "transfer load in step " << Step() << IO::endl;
}

bool PARTICLEALGORITHM::ParticleAlgorithm::CheckLoadRedistributionNeeded()
{
  bool redistributeload = writerestartthisstep_;

  // percentage limit
  const double percentagelimit = 0.1;

  // get number of particles on this processor
  int numberofparticles = particleengine_->GetNumberOfParticles();

  // percentage change of particles on this processor
  double percentagechange = 0.0;
  if (numparticlesafterlastloadbalance_ > 0)
    percentagechange =
        std::abs(static_cast<double>(numberofparticles - numparticlesafterlastloadbalance_) /
                 numparticlesafterlastloadbalance_);

  // get maximum percentage change of particles
  double maxpercentagechange = 0.0;
  Comm().MaxAll(&percentagechange, &maxpercentagechange, 1);

  // criterion for load redistribution based on maximum percentage change of the number of particles
  redistributeload |= (maxpercentagechange > percentagelimit);

  return redistributeload;
}

void PARTICLEALGORITHM::ParticleAlgorithm::DistributeLoadAmongProcs()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::ParticleAlgorithm::DistributeLoadAmongProcs");

  // dynamic load balancing
  particleengine_->DynamicLoadBalancing();

  // get number of particles on this processor
  numparticlesafterlastloadbalance_ = particleengine_->GetNumberOfParticles();

  if (particlewall_)
  {
    // update bin row and column map
    particlewall_->UpdateBinRowAndColMap(
        particleengine_->GetBinRowMap(), particleengine_->GetBinColMap());

    // distribute wall elements and nodes
    particlewall_->DistributeWallElementsAndNodes();
  }

  // communicate rigid body
  if (particlerigidbody_) particlerigidbody_->CommunicateRigidBody();

  // communicate interaction history
  if (particleinteraction_) particleinteraction_->CommunicateInteractionHistory();

  // short screen output
  if (myrank_ == 0) IO::cout(IO::verbose) << "distribute load in step " << Step() << IO::endl;
}

void PARTICLEALGORITHM::ParticleAlgorithm::BuildPotentialNeighborRelation()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::ParticleAlgorithm::BuildPotentialNeighborRelation");

  // build particle to particle neighbors
  particleengine_->BuildParticleToParticleNeighbors();

  if (particlewall_)
  {
    // relate bins to column wall elements
    particlewall_->RelateBinsToColWallEles();

    // build particle to wall neighbors
    particlewall_->BuildParticleToWallNeighbors(particleengine_->GetParticlesToBins());
  }
}

void PARTICLEALGORITHM::ParticleAlgorithm::SetInitialConditions()
{
  // create and init particle initial field handler
  std::unique_ptr<PARTICLEALGORITHM::InitialFieldHandler> initialfield =
      std::unique_ptr<PARTICLEALGORITHM::InitialFieldHandler>(
          new PARTICLEALGORITHM::InitialFieldHandler(params_));
  initialfield->Init();

  // setup particle initial field handler
  initialfield->Setup(particleengine_);

  // set initial fields
  initialfield->SetInitialFields();
}

void PARTICLEALGORITHM::ParticleAlgorithm::SetCurrentTime()
{
  // set current time in particle time integration
  particletimint_->SetCurrentTime(Time());

  // set current time in particle interaction
  if (particleinteraction_) particleinteraction_->SetCurrentTime(Time());
}

void PARTICLEALGORITHM::ParticleAlgorithm::SetCurrentStepSize()
{
  // set current step size in particle interaction
  if (particleinteraction_) particleinteraction_->SetCurrentStepSize(Dt());
}

void PARTICLEALGORITHM::ParticleAlgorithm::SetCurrentWriteResultFlag()
{
  // set current write result flag in particle interaction
  if (particleinteraction_) particleinteraction_->SetCurrentWriteResultFlag(writeresultsthisstep_);
}

void PARTICLEALGORITHM::ParticleAlgorithm::EvaluateTimeStep()
{
  // clear forces and torques
  if (particlerigidbody_) particlerigidbody_->ClearForcesAndTorques();

  // set gravity acceleration
  if (particlegravity_) SetGravityAcceleration();

  // evaluate particle interactions
  if (particleinteraction_) particleinteraction_->EvaluateInteractions();

  // apply viscous damping contribution
  if (viscousdamping_) viscousdamping_->ApplyViscousDamping();

  // compute accelerations of rigid bodies
  if (particlerigidbody_ and particleinteraction_) particlerigidbody_->ComputeAccelerations();
}

void PARTICLEALGORITHM::ParticleAlgorithm::SetGravityAcceleration()
{
  std::vector<double> scaled_gravity(3);

  // get gravity acceleration
  particlegravity_->GetGravityAcceleration(Time(), scaled_gravity);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengine_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle->GetParticleTypes())
  {
    // gravity is not set for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // gravity is not set for open boundary particles
    if (typeEnum == PARTICLEENGINE::DirichletPhase or typeEnum == PARTICLEENGINE::NeumannPhase)
      continue;

    // set gravity acceleration for all particles of current type
    particlecontainerbundle->SetStateSpecificContainer(
        scaled_gravity, PARTICLEENGINE::Acceleration, typeEnum);
  }

  // add gravity acceleration
  if (particlerigidbody_) particlerigidbody_->AddGravityAcceleration(scaled_gravity);

  // set scaled gravity in particle interaction handler
  if (particleinteraction_) particleinteraction_->SetGravity(scaled_gravity);
}

BACI_NAMESPACE_CLOSE
