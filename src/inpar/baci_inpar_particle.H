/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle problems

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef BACI_INPAR_PARTICLE_H
#define BACI_INPAR_PARTICLE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "baci_inpar_parameterlist_utils.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*---------------------------------------------------------------------------*
 | input parameters for particle problems                     sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace INPAR
{
  namespace PARTICLE
  {
    /*---------------------------------------------------------------------------*
     | general control parameters for particle simulations                       |
     *---------------------------------------------------------------------------*/

    //! type of particle time integration
    enum DynamicType
    {
      dyna_semiimpliciteuler,  //! semi implicit euler scheme (explicit)
      dyna_velocityverlet      //! velocity verlet scheme (explicit)
    };

    //! type of particle interaction
    enum InteractionType
    {
      interaction_none,  //! no particle interaction
      interaction_sph,   //! smoothed particle hydrodynamics
      interaction_dem    //! discrete element method
    };

    //! data format for written numeric data via vtp
    enum OutputDataFormat
    {
      binary,
      ascii
    };

    //! type of particle wall source
    enum ParticleWallSource
    {
      NoParticleWall,    //! no particle wall
      DiscretCondition,  //! particle wall from discretization condition
      BoundingBox        //! particle wall from bounding box
    };

    /*---------------------------------------------------------------------------*
     | smoothed particle hydrodynamics (SPH) specific control parameters         |
     *---------------------------------------------------------------------------*/

    //! type of smoothed particle hydrodynamics kernel
    enum KernelType
    {
      CubicSpline,
      QuinticSpline
    };

    //! kernel space dimension number
    enum KernelSpaceDimension
    {
      Kernel1D,
      Kernel2D,
      Kernel3D
    };

    //! type of smoothed particle hydrodynamics equation of state
    enum EquationOfStateType
    {
      GenTait,
      IdealGas
    };

    //! type of smoothed particle hydrodynamics momentum formulation
    enum MomentumFormulationType
    {
      AdamiMomentumFormulation,
      MonaghanMomentumFormulation
    };

    //! type of density evaluation scheme
    enum DensityEvaluationScheme
    {
      DensitySummation,
      DensityIntegration,
      DensityPredictCorrect
    };

    //! type of density correction scheme
    enum DensityCorrectionScheme
    {
      NoCorrection,
      InteriorCorrection,
      NormalizedCorrection,
      RandlesCorrection
    };

    //! type of boundary particle formulation
    enum BoundaryParticleFormulationType
    {
      NoBoundaryFormulation,
      AdamiBoundaryFormulation
    };

    //! type of boundary particle interaction
    enum BoundaryParticleInteraction
    {
      NoSlipBoundaryParticle,
      FreeSlipBoundaryParticle
    };

    //! type of wall formulation
    enum WallFormulationType
    {
      NoWallFormulation,
      VirtualParticleWallFormulation
    };

    //! type of transport velocity formulation
    enum TransportVelocityFormulation
    {
      NoTransportVelocity,
      StandardTransportVelocity,
      GeneralizedTransportVelocity
    };

    //! type of temperature evaluation scheme
    enum TemperatureEvaluationScheme
    {
      NoTemperatureEvaluation,
      TemperatureIntegration
    };

    //! type of heat source
    enum HeatSourceType
    {
      NoHeatSource,
      VolumeHeatSource,
      SurfaceHeatSource
    };

    //! type of surface tension formulation
    enum SurfaceTensionFormulation
    {
      NoSurfaceTension,
      ContinuumSurfaceForce
    };

    //! type of dirichlet open boundary
    enum DirichletOpenBoundaryType
    {
      NoDirichletOpenBoundary,
      DirichletNormalToPlane
    };

    //! type of neumann open boundary
    enum NeumannOpenBoundaryType
    {
      NoNeumannOpenBoundary,
      NeumannNormalToPlane
    };

    //! type of phase change
    enum PhaseChangeType
    {
      NoPhaseChange,
      OneWayScalarBelowToAbovePhaseChange,
      OneWayScalarAboveToBelowPhaseChange,
      TwoWayScalarPhaseChange
    };

    //! type of rigid particle contact
    enum RigidParticleContactType
    {
      NoRigidParticleContact,
      ElasticRigidParticleContact
    };

    /*---------------------------------------------------------------------------*
     | discrete element method (DEM) specific control parameters                 |
     *---------------------------------------------------------------------------*/

    //! type of normal contact law
    enum NormalContact
    {
      NormalLinSpring,
      NormalLinSpringDamp,
      NormalHertz,
      NormalLeeHerrmann,
      NormalKuwabaraKono,
      NormalTsuji,
      NormalForceDistanceCurve,
      NormalForceDistanceCurveDamp
    };

    //! type of tangential contact law
    enum TangentialContact
    {
      NoTangentialContact,
      TangentialLinSpringDamp
    };

    //! type of rolling contact law
    enum RollingContact
    {
      NoRollingContact,
      RollingViscous,
      RollingCoulomb
    };

    //! type of adhesion law
    enum AdhesionLaw
    {
      NoAdhesion,
      AdhesionVdWDMT,
      AdhesionRegDMT
    };

    //! type of (random) surface energy distribution
    enum SurfaceEnergyDistribution
    {
      ConstantSurfaceEnergy,
      NormalSurfaceEnergyDistribution,
      LogNormalSurfaceEnergyDistribution
    };

    //! type of initial particle radius assignment
    enum InitialRadiusAssignment
    {
      RadiusFromParticleMaterial,
      RadiusFromParticleInput,
      NormalRadiusDistribution,
      LogNormalRadiusDistribution
    };

    //! set the particle parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    //! set the particle conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace PARTICLE

}  // namespace INPAR

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
