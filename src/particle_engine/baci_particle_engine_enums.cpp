/*---------------------------------------------------------------------------*/
/*! \file
\brief defining enums for particle problem
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_engine_enums.H"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
int PARTICLEENGINE::EnumToStateDim(const enum ParticleState& state)
{
  int dim = 0;

  switch (state)
  {
    // scalar states
    case Radius:
    case InitRadius:
    case Mass:
    case Density:
    case Pressure:
    case Temperature:
    case RigidBodyColor:
    case Inertia:
    case DensitySum:
    case DensityDot:
    case TemperatureDot:
    case BoundaryPressure:
    case Colorfield:
    case Curvature:
    case WallColorfield:
    case LastIterDensity:
    case LastIterTemperature:
      dim = 1;
      break;

    // vectorial states
    case Position:
    case Velocity:
    case Acceleration:
    case LastTransferPosition:
    case ModifiedVelocity:
    case ModifiedAcceleration:
    case ReferencePosition:
    case RelativePosition:
    case RelativePositionBodyFrame:
    case BoundaryVelocity:
    case ColorfieldGradient:
    case InterfaceNormal:
    case WallInterfaceNormal:
    case TemperatureGradient:
    case AngularVelocity:
    case AngularAcceleration:
    case Force:
    case Moment:
    case LastIterPosition:
    case LastIterVelocity:
    case LastIterAcceleration:
    case LastIterAngularVelocity:
    case LastIterAngularAcceleration:
    case LastIterModifiedAcceleration:
      dim = 3;
      break;

    default:
      dserror("particle state unknown!");
  }

  return dim;
}

std::string PARTICLEENGINE::EnumToStateName(const enum ParticleState& state)
{
  std::string name;

  switch (state)
  {
    case Radius:
      name = "radius";
      break;
    case InitRadius:
      name = "initial radius";
      break;
    case Mass:
      name = "mass";
      break;
    case Density:
      name = "density";
      break;
    case DensitySum:
      name = "density sum";
      break;
    case DensityDot:
      name = "density dot";
      break;
    case Pressure:
      name = "pressure";
      break;
    case Temperature:
      name = "temperature";
      break;
    case RigidBodyColor:
      name = "rigid body color";
      break;
    case Inertia:
      name = "inertia";
      break;
    case TemperatureDot:
      name = "temperature dot";
      break;
    case Position:
      name = "position";
      break;
    case Velocity:
      name = "velocity";
      break;
    case Acceleration:
      name = "acceleration";
      break;
    case AngularVelocity:
      name = "angular velocity";
      break;
    case AngularAcceleration:
      name = "angular acceleration";
      break;
    case Force:
      name = "force";
      break;
    case Moment:
      name = "moment";
      break;
    case LastTransferPosition:
      name = "position last transfer";
      break;
    case ReferencePosition:
      name = "reference position";
      break;
    case RelativePosition:
      name = "relative position";
      break;
    case RelativePositionBodyFrame:
      name = "body frame relative position";
      break;
    case ModifiedVelocity:
      name = "modified velocity";
      break;
    case ModifiedAcceleration:
      name = "modified acceleration";
      break;
    case BoundaryPressure:
      name = "boundary pressure";
      break;
    case BoundaryVelocity:
      name = "boundary velocity";
      break;
    case Colorfield:
      name = "colorfield";
      break;
    case ColorfieldGradient:
      name = "colorfield gradient";
      break;
    case InterfaceNormal:
      name = "interface normal";
      break;
    case Curvature:
      name = "curvature";
      break;
    case WallColorfield:
      name = "wall colorfield";
      break;
    case WallInterfaceNormal:
      name = "wall interface normal";
      break;
    case TemperatureGradient:
      name = "temperature gradient";
      break;
    case LastIterPosition:
      name = "position last iteration";
      break;
    case LastIterVelocity:
      name = "velocity last iteration";
      break;
    case LastIterAcceleration:
      name = "acceleration last iteration";
      break;
    case LastIterAngularVelocity:
      name = "angular velocity last iteration";
      break;
    case LastIterAngularAcceleration:
      name = "angular acceleration last iteration";
      break;
    case LastIterModifiedAcceleration:
      name = "modified acceleration last iteration";
      break;
    case LastIterDensity:
      name = "density last iteration";
      break;
    case LastIterTemperature:
      name = "temperature last iteration";
      break;
    default:
      dserror("particle state unknown!");
  }

  return name;
}

enum PARTICLEENGINE::ParticleState PARTICLEENGINE::EnumFromStateName(const std::string& name)
{
  enum ParticleState state;

  if (name == "density")
    state = Density;
  else if (name == "pressure")
    state = Pressure;
  else if (name == "temperature")
    state = Temperature;
  else
    dserror("particle state '%s' unknown!", name.c_str());

  return state;
}

std::string PARTICLEENGINE::EnumToTypeName(const enum ParticleType& type)
{
  std::string name;

  switch (type)
  {
    case Phase1:
      name = "phase1";
      break;
    case Phase2:
      name = "phase2";
      break;
    case BoundaryPhase:
      name = "boundaryphase";
      break;
    case RigidPhase:
      name = "rigidphase";
      break;
    case DirichletPhase:
      name = "dirichletphase";
      break;
    case NeumannPhase:
      name = "neumannphase";
      break;
    default:
      dserror("particle type unknown!");
  }

  return name;
}

enum PARTICLEENGINE::ParticleType PARTICLEENGINE::EnumFromTypeName(const std::string& name)
{
  enum ParticleType type;

  if (name == "phase1")
    type = Phase1;
  else if (name == "phase2")
    type = Phase2;
  else if (name == "boundaryphase")
    type = BoundaryPhase;
  else if (name == "rigidphase")
    type = RigidPhase;
  else if (name == "dirichletphase")
    type = DirichletPhase;
  else if (name == "neumannphase")
    type = NeumannPhase;
  else
    dserror("particle type '%s' unknown!", name.c_str());

  return type;
}

std::string PARTICLEENGINE::EnumToStatusName(const enum ParticleStatus& status)
{
  std::string name;

  switch (status)
  {
    case Owned:
      name = "owned";
      break;
    case Ghosted:
      name = "ghosted";
      break;
    default:
      dserror("particle status unknown!");
  }

  return name;
}

BACI_NAMESPACE_CLOSE
