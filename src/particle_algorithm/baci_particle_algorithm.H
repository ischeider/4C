/*---------------------------------------------------------------------------*/
/*! \file
\brief algorithm to control particle simulations
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_ALGORITHM_H
#define BACI_PARTICLE_ALGORITHM_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_adapter_algorithmbase.H"
#include "baci_particle_engine.H"
#include "baci_particle_engine_typedefs.H"
#include "baci_particle_wall.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class TimInt;
  class GravityHandler;
  class ViscousDampingHandler;
  class ParticleGrowthHandler;
}  // namespace PARTICLEALGORITHM

namespace PARTICLEINTERACTION
{
  class ParticleInteractionBase;
}

namespace PARTICLEENGINE
{
  class ParticleObject;
}

namespace PARTICLERIGIDBODY
{
  class RigidBodyHandler;
}

namespace DRT
{
  class ResultTest;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief algorithm to control particle simulations
   *
   * The particle algorithm controls all major steps of particle simulations. Besides initialization
   * and setup of the problem this mainly concerns the control of the time loop (with the
   * integration of the particle step) consisting of the following steps:
   * - update and maintain valid connectivity (parallel distribution, load transfer, ghosting,
   *   neighbor pair relations)
   * - evaluation of particle accelerations from gravity acceleration, interactions, and viscous
   *   damping
   * - (runtime) output of fields
   * - write restart information
   *
   * \author Sebastian Fuchs \date 04/2018
   */

  class ParticleAlgorithm final : public ADAPTER::AlgorithmBase
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] comm   communicator
     * \param[in] params particle simulation parameter list
     */
    explicit ParticleAlgorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~ParticleAlgorithm() override;

    /*!
     * \brief init particle algorithm
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] initialparticles particle objects read in from restart
     */
    void Init(std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles);

    /*!
     * \brief setup particle algorithm
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void Setup();

    /*!
     * \brief read restart information for given time step
     *
     * Read the restart information in the same order as it is written.
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] restartstep restart step
     */
    void ReadRestart(const int restartstep) override;

    /*!
     * \brief time loop for particle problem
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void Timeloop();

    /*!
     * \brief prepare time step
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] print_header flag to control output of time step header
     */
    void PrepareTimeStep(bool print_header = true);

    /*!
     * \brief pre evaluate time step
     *
     * \author Sebastian Fuchs \date 11/2020
     */
    void PreEvaluateTimeStep();

    /*!
     * \brief integrate time step
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void IntegrateTimeStep();

    /*!
     * \brief post evaluate time step
     *
     * \author Sebastian Fuchs \date 11/2019
     */
    void PostEvaluateTimeStep();

    /*!
     * \brief write output
     *
     * \author Sebastian Fuchs \date 07/2020
     */
    void WriteOutput() const;

    /*!
     * \brief write restart information
     *
     * The order of writing the restart information in the particle framework is important: First
     * restart information is written by the bin discretization writer, afterwards restart
     * information is written by the wall discretization writer. This order is important, otherwise
     * the restart fails.
     *
     * \author Sebastian Fuchs \date 07/2020
     */
    void WriteRestart() const;

    /*!
     * \brief create particle field specific result test objects
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \return particle field specific result test objects
     */
    std::vector<std::shared_ptr<DRT::ResultTest>> CreateResultTests();

    /*!
     * \brief get interface to particle engine
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \return interface to particle engine
     */
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> GetParticleEngineInterface() const
    {
      return particleengine_;
    }

    /*!
     * \brief get interface to particle wall handler
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \return interface to particle wall handler
     */
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> GetParticleWallHandlerInterface() const
    {
      return particlewall_;
    }

   private:
    //! \name init and setup methods
    //! @{

    /*!
     * \brief init particle engine
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void InitParticleEngine();

    /*!
     * \brief init particle wall handler
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void InitParticleWall();

    /*!
     * \brief init rigid body handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void InitParticleRigidBody();

    /*!
     * \brief init particle time integration
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void InitParticleTimeIntegration();

    /*!
     * \brief init particle interaction handler
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void InitParticleInteraction();

    /*!
     * \brief init particle gravity handler
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void InitParticleGravity();

    /*!
     * \brief init particle growth handler
     *
     * \author Ingo Scheider \date 10/2023
     */
    void InitParticleGrowth();

    /*!
     * \brief init viscous damping handler
     *
     * \author Sebastian Fuchs \date 02/2019
     */
    void InitViscousDamping();

    /*!
     * \brief generate initial particles
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void GenerateInitialParticles();

    /*!
     * \brief determine all particle types
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void DetermineParticleTypes();

    /*!
     * \brief determine particle states of all particle types
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void DetermineParticleStatesOfParticleTypes();

    /*!
     * \brief setup initial particles
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void SetupInitialParticles();

    /*!
     * \brief setup initial rigid bodies
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void SetupInitialRigidBodies();

    /*!
     * \brief setup initial states
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void SetupInitialStates();

    //! @}

    /*!
     * \brief update and maintain valid connectivity in particle problems
     *
     * The main control routine of particle problems to maintain valid particle connectivity:
     * - parallel distribution of particles
     * - transfer of particles to new bins
     * - ghosting of particles
     * - neighbor pair relations
     *
     * This method includes several safety checks that are part of the called sub-methods.
     *
     * \author Sebastian Fuchs \date 06/2018
     */
    void UpdateConnectivity();

    /*!
     * \brief check load transfer
     *
     * Check if a load transfer is needed based on the following criteria:
     * - maximum position increment of particles or walls
     * - invalid particle connectivity
     * - invalid particle neighbor pair relation
     * - invalid particle wall neighbor pair relation
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \return flag indicating load transfer is needed
     */
    bool CheckLoadTransferNeeded();

    /*!
     * \brief check maximum position increment
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \return flag indicating load transfer is needed due to maximum position increment
     */
    bool CheckMaxPositionIncrement();

    /*!
     * \brief get maximum particle position increment since last transfer
     *
     * \author Sebastian Fuchs \date 06/2018
     *
     * \param allprocmaxpositionincrement maximum particle position increment
     */
    void GetMaxParticlePositionIncrement(double& allprocmaxpositionincrement);

    /*!
     * \brief transfer load between processors
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    void TransferLoadBetweenProcs();

    /*!
     * \brief check load redistribution
     *
     * Check if a load redistribution is needed based on the difference in the absolute number of
     * particles since the last redistribution.
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \return flag indicating load redistribution is needed
     */
    bool CheckLoadRedistributionNeeded();

    /*!
     * \brief distribute load among processors
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    void DistributeLoadAmongProcs();

    /*!
     * \brief build potential neighbor relation
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    void BuildPotentialNeighborRelation();

    /*!
     * \brief set initial conditions
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void SetInitialConditions();

    /*!
     * \brief set current time
     *
     * \author Sebastian Fuchs \date 08/2018
     */
    void SetCurrentTime();

    /*!
     * \brief set current step size
     *
     * \author Sebastian Fuchs \date 08/2018
     */
    void SetCurrentStepSize();

    /*!
     * \brief set current write result flag
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    void SetCurrentWriteResultFlag();

    /*!
     * \brief evaluate time step
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void EvaluateTimeStep();

    /*!
     * \brief set gravity acceleration
     *
     * \author Sebastian Fuchs \date 06/2018
     */
    void SetGravityAcceleration();

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngine> particleengine_;

    //! particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerBase> particlewall_;

    //! rigid body handler
    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyHandler> particlerigidbody_;

    //! particle time integration
    std::unique_ptr<PARTICLEALGORITHM::TimInt> particletimint_;

    //! particle interaction
    std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionBase> particleinteraction_;

    //! particle growth handler
    std::unique_ptr<PARTICLEALGORITHM::ParticleGrowthHandler> particlegrowth_;

    //! particle gravity handler
    std::unique_ptr<PARTICLEALGORITHM::GravityHandler> particlegravity_;

    //! particle viscous damping handler
    std::unique_ptr<PARTICLEALGORITHM::ViscousDampingHandler> viscousdamping_;

    //! map of particle types and corresponding states
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes_;

    //! vector of initial or generated particles to distribute
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr> particlestodistribute_;

    //! number of particles on this processor after last load balance
    int numparticlesafterlastloadbalance_;

    //! transfer particles to new bins every time step
    bool transferevery_;

    //! write results interval
    const int writeresultsevery_;

    //! write restart interval
    const int writerestartevery_;

    //! result control flag
    bool writeresultsthisstep_;

    //! restart control flag
    bool writerestartthisstep_;

    //! simulation is restarted
    bool isrestarted_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
