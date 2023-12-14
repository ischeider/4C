/*---------------------------------------------------------------------------*/
/*! \file
\brief particle growth handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_ALGORITHM_GROWTH_H
#define BACI_PARTICLE_ALGORITHM_GROWTH_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_engine_typedefs.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief particle growth handler for particle simulations
   */
  class ParticleGrowthHandler
  {
   public:
    /*!
     * \param[in] growthfunc function for particle growth
     */
    explicit ParticleGrowthHandler(const int growthfunc);

    /*!
     * \brief init particle growth handler
     */
    void Init();

    /*!
     * \brief setup particle growth handler
     */
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);


    //! insert state (InitialRadius) of particle type with growth
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes);

    //! set initial states
    void SetInitialStates();

    /*!
     * \brief apply viscous damping contribution
     */
    void PreEvaluateTimeStep(double currenttime);

   private:
    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle growth function
    const int growthfunc_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif

BACI_NAMESPACE_CLOSE