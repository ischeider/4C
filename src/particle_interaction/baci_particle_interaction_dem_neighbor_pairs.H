/*---------------------------------------------------------------------------*/
/*! \file
\brief neighbor pair handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_DEM_NEIGHBOR_PAIRS_H
#define BACI_PARTICLE_INTERACTION_DEM_NEIGHBOR_PAIRS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_engine_enums.H"
#include "baci_particle_engine_typedefs.H"
#include "baci_particle_interaction_dem_neighbor_pair_struct.H"

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

/*---------------------------------------------------------------------------*
 | type definitions                                                          |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  using DEMParticlePairData = std::vector<PARTICLEINTERACTION::DEMParticlePair>;
  using DEMParticleWallPairData = std::vector<PARTICLEINTERACTION::DEMParticleWallPair>;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMNeighborPairs final
  {
   public:
    //! constructor
    explicit DEMNeighborPairs();

    //! init neighbor pair handler
    void Init();

    //! setup neighbor pair handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface);

    //! get reference to particle pair data
    inline const DEMParticlePairData& GetRefToParticlePairData() const
    {
      return particlepairdata_;
    };

    //! get reference to particle-wall pair data
    inline const DEMParticleWallPairData& GetRefToParticleWallPairData() const
    {
      return particlewallpairdata_;
    };

    //! get reference to adhesion particle pair data
    inline const DEMParticlePairData& GetRefToParticlePairAdhesionData() const
    {
      return particlepairadhesiondata_;
    };

    //! get reference to adhesion particle-wall pair data
    inline const DEMParticleWallPairData& GetRefToParticleWallPairAdhesionData() const
    {
      return particlewallpairadhesiondata_;
    };

    //! evaluate neighbor pairs
    void EvaluateNeighborPairs(const double& interaction_horizon);

    //! evaluate adhesion neighbor pairs
    void EvaluateNeighborPairsAdhesion(const double& adhesion_distance);

   private:
    //! evaluate particle pairs
    void EvaluateParticlePairs(const double& interaction_horizon);

    //! evaluate particle-wall pairs
    void EvaluateParticleWallPairs(const double& interaction_horizon);

    //! evaluate adhesion particle pairs
    void EvaluateParticlePairsAdhesion(const double& adhesion_distance);

    //! evaluate adhesion particle-wall pairs
    void EvaluateParticleWallPairsAdhesion(const double& adhesion_distance);

    //! particle pair data with evaluated quantities
    DEMParticlePairData particlepairdata_;

    //! particle-wall pair data with evaluated quantities
    DEMParticleWallPairData particlewallpairdata_;

    //! adhesion particle pair data with evaluated quantities
    DEMParticlePairData particlepairadhesiondata_;

    //! adhesion particle-wall pair data with evaluated quantities
    DEMParticleWallPairData particlewallpairadhesiondata_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
