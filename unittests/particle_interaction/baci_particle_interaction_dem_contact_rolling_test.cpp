/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for rolling contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_particle_interaction_dem_contact_rolling.H"

#include "baci_inpar_validparameters.H"
#include "baci_particle_interaction_utils.H"
#include "baci_unittest_utils_assertions_test.H"

namespace
{
  using namespace BACI;

  class DEMContactRollingViscousTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactRollingViscous> contactrolling_;

    const double e_ = 0.8;
    const double nue_ = 0.4;
    const double mu_rolling_ = 0.2;

    const double young_ = 200.0e3;
    const double v_max_ = 0.025;

    const double k_normal_per_radius_ = 4.0;

    DEMContactRollingViscousTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("COEFF_RESTITUTION", e_);
      params_dem.set("POISSON_RATIO", nue_);
      params_dem.set("FRICT_COEFF_ROLL", mu_rolling_);

      params_dem.set("YOUNG_MODULUS", young_);
      params_dem.set("MAX_VELOCITY", v_max_);

      // create rolling contact handler
      contactrolling_ = std::make_unique<PARTICLEINTERACTION::DEMContactRollingViscous>(params_dem);

      // init rolling contact handler
      contactrolling_->Init();

      // setup rolling contact handler
      contactrolling_->Setup(k_normal_per_radius_);
    }
    // note: the public functions Init() and Setup() of class DEMContactRollingViscous are
    // called in the constructor and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactRollingViscousTest, EffectiveRadiusParticle)
  {
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double gap = -0.3;

    double r_eff = 0.0;
    contactrolling_->EffectiveRadiusParticle(&rad_i, &rad_j, gap, r_eff);

    const double r_eff_ref = rad_i * rad_j / (rad_i + rad_j);

    EXPECT_NEAR(r_eff, r_eff_ref, 1.0e-12);
  }

  TEST_F(DEMContactRollingViscousTest, EffectiveRadiusParticleNullptrRadJ)
  {
    const double rad_i = 1.2;
    const double gap = -0.3;

    double r_eff = 0.0;
    contactrolling_->EffectiveRadiusParticle(&rad_i, nullptr, gap, r_eff);

    const double r_eff_ref = rad_i;

    EXPECT_NEAR(r_eff, r_eff_ref, 1.0e-12);
  }

  TEST_F(DEMContactRollingViscousTest, RelativeRollingVelocity)
  {
    const double r_eff = 0.5;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double angvel_i[3] = {0.0};
    angvel_i[0] = -0.03;
    angvel_i[1] = 0.1;
    angvel_i[2] = 0.12;

    double angvel_j[3] = {0.0};
    angvel_j[0] = 0.15;
    angvel_j[1] = -0.2;
    angvel_j[2] = 0.0;

    double v_rel_rolling[3] = {0.0};
    contactrolling_->RelativeRollingVelocity(r_eff, e_ji, angvel_i, angvel_j, v_rel_rolling);

    double v_rel_rolling_ref[3] = {0.0};
    PARTICLEINTERACTION::UTILS::VecSetCross(v_rel_rolling_ref, angvel_i, e_ji);
    PARTICLEINTERACTION::UTILS::VecAddCross(v_rel_rolling_ref, e_ji, angvel_j);

    BACI_EXPECT_ITERABLE_NEAR(v_rel_rolling, v_rel_rolling_ref, 3, 1.0e-12);
  }

  TEST_F(DEMContactRollingViscousTest, RelativeRollingVelocityNullptrAngvelJ)
  {
    const double r_eff = 0.5;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double angvel_i[3] = {0.0};
    angvel_i[0] = -0.03;
    angvel_i[1] = 0.1;
    angvel_i[2] = 0.12;

    double v_rel_rolling[3] = {0.0};
    contactrolling_->RelativeRollingVelocity(r_eff, e_ji, angvel_i, nullptr, v_rel_rolling);

    double v_rel_rolling_ref[3] = {0.0};
    PARTICLEINTERACTION::UTILS::VecSetCross(v_rel_rolling_ref, angvel_i, e_ji);

    BACI_EXPECT_ITERABLE_NEAR(v_rel_rolling, v_rel_rolling_ref, 3, 1.0e-12);
  }

  TEST_F(DEMContactRollingViscousTest, RollingContactMoment)
  {
    double gap_rolling[3] = {0.0};
    gap_rolling[0] = 0.1;
    gap_rolling[1] = 0.05;
    gap_rolling[2] = -0.25;

    bool stick_rolling = false;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double v_rel_rolling[3] = {0.0};
    v_rel_rolling[0] = -0.03;
    v_rel_rolling[1] = 0.1;
    v_rel_rolling[2] = 0.12;

    const double m_eff = 2.5;
    const double r_eff = 0.5;
    const double normalcontactforce = 1.5;

    double rollingcontactmoment[3] = {0.0};
    contactrolling_->RollingContactMoment(gap_rolling, stick_rolling, e_ji, v_rel_rolling, m_eff,
        r_eff, mu_rolling_, normalcontactforce, rollingcontactmoment);

    double rollingcontactmoment_ref[3] = {0.0};
    const double fac = young_ / (1.0 - PARTICLEINTERACTION::UTILS::Pow<2>(nue_));
    const double c_1 = 1.15344;
    const double d_rolling_fac =
        mu_rolling_ * (1.0 - e_) / (c_1 * std::pow(fac, 0.4) * std::pow(v_max_, 0.2));
    const double d_rolling = d_rolling_fac * std::pow(0.5 * r_eff, -0.2);

    double rollingcontactforce[3];
    PARTICLEINTERACTION::UTILS::VecSetScale(
        rollingcontactforce, -(d_rolling * normalcontactforce), v_rel_rolling);

    PARTICLEINTERACTION::UTILS::VecSetCross(rollingcontactmoment_ref, rollingcontactforce, e_ji);
    PARTICLEINTERACTION::UTILS::VecScale(rollingcontactmoment_ref, r_eff);

    for (int i = 0; i < 3; ++i)
      EXPECT_NEAR(rollingcontactmoment[i], rollingcontactmoment_ref[i], 1.0e-12);
  }

  TEST_F(DEMContactRollingViscousTest, RollingPotentialEnergy)
  {
    double gap_rolling[3] = {0.0};

    const double rollingpotentialenergy_ref = 0.0;

    double rollingpotentialenergy = 0.0;
    contactrolling_->RollingPotentialEnergy(gap_rolling, rollingpotentialenergy);

    EXPECT_NEAR(rollingpotentialenergy, rollingpotentialenergy_ref, 1.0e-12);
  }

  class DEMContactRollingCoulombTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactRollingCoulomb> contactrolling_;

    const double e_ = 0.8;
    const double nue_ = 0.4;
    const double mu_rolling_ = 0.2;

    const double k_normal_per_radius_ = 4.0;

    DEMContactRollingCoulombTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("COEFF_RESTITUTION", e_);
      params_dem.set("POISSON_RATIO", nue_);
      params_dem.set("FRICT_COEFF_ROLL", mu_rolling_);

      // create rolling contact handler
      contactrolling_ = std::make_unique<PARTICLEINTERACTION::DEMContactRollingCoulomb>(params_dem);

      // init rolling contact handler
      contactrolling_->Init();

      // setup rolling contact handler
      contactrolling_->Setup(k_normal_per_radius_);
    }
    // note: the public functions Init() and Setup() of class DEMContactRollingViscous are
    // called in the constructor and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactRollingCoulombTest, EffectiveRadiusParticle)
  {
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double gap = -0.3;

    double r_eff = 0.0;
    contactrolling_->EffectiveRadiusParticle(&rad_i, &rad_j, gap, r_eff);

    const double r_eff_ref = (rad_i + 0.5 * gap) * (rad_j + 0.5 * gap) / (rad_i + rad_j + gap);

    EXPECT_NEAR(r_eff, r_eff_ref, 1.0e-12);
  }

  TEST_F(DEMContactRollingCoulombTest, EffectiveRadiusParticleNullptrRadJ)
  {
    const double rad_i = 1.2;
    const double gap = -0.3;

    double r_eff = 0.0;
    contactrolling_->EffectiveRadiusParticle(&rad_i, nullptr, gap, r_eff);

    const double r_eff_ref = rad_i + gap;

    EXPECT_NEAR(r_eff, r_eff_ref, 1.0e-12);
  }

  TEST_F(DEMContactRollingCoulombTest, RelativeRollingVelocity)
  {
    const double r_eff = 0.5;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double angvel_i[3] = {0.0};
    angvel_i[0] = -0.03;
    angvel_i[1] = 0.1;
    angvel_i[2] = 0.12;

    double angvel_j[3] = {0.0};
    angvel_j[0] = 0.15;
    angvel_j[1] = -0.2;
    angvel_j[2] = 0.0;

    double v_rel_rolling[3] = {0.0};
    contactrolling_->RelativeRollingVelocity(r_eff, e_ji, angvel_i, angvel_j, v_rel_rolling);

    double v_rel_rolling_ref[3] = {0.0};
    PARTICLEINTERACTION::UTILS::VecSetCross(v_rel_rolling_ref, e_ji, angvel_i);
    PARTICLEINTERACTION::UTILS::VecAddCross(v_rel_rolling_ref, angvel_j, e_ji);
    PARTICLEINTERACTION::UTILS::VecScale(v_rel_rolling_ref, r_eff);

    BACI_EXPECT_ITERABLE_NEAR(v_rel_rolling, v_rel_rolling_ref, 3, 1.0e-12);
  }

  TEST_F(DEMContactRollingCoulombTest, RelativeRollingVelocityNullptrAngvelJ)
  {
    const double r_eff = 0.5;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double angvel_i[3] = {0.0};
    angvel_i[0] = -0.03;
    angvel_i[1] = 0.1;
    angvel_i[2] = 0.12;

    double v_rel_rolling[3] = {0.0};
    contactrolling_->RelativeRollingVelocity(r_eff, e_ji, angvel_i, nullptr, v_rel_rolling);

    double v_rel_rolling_ref[3] = {0.0};
    PARTICLEINTERACTION::UTILS::VecSetCross(v_rel_rolling_ref, e_ji, angvel_i);
    PARTICLEINTERACTION::UTILS::VecScale(v_rel_rolling_ref, r_eff);

    BACI_EXPECT_ITERABLE_NEAR(v_rel_rolling, v_rel_rolling_ref, 3, 1.0e-12);
  }

  TEST_F(DEMContactRollingCoulombTest, RollingContactMomentStick)
  {
    double gap_rolling[3] = {0.0};
    gap_rolling[0] = 0.1;
    gap_rolling[1] = 0.05;
    gap_rolling[2] = -0.25;

    bool stick_rolling = true;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double v_rel_rolling[3] = {0.0};
    v_rel_rolling[0] = -0.03;
    v_rel_rolling[1] = 0.1;
    v_rel_rolling[2] = 0.12;

    const double m_eff = 2.5;
    const double r_eff = 0.5;
    const double normalcontactforce = 1.5;

    double rollingcontactmoment[3] = {0.0};
    contactrolling_->RollingContactMoment(gap_rolling, stick_rolling, e_ji, v_rel_rolling, m_eff,
        r_eff, mu_rolling_, normalcontactforce, rollingcontactmoment);

    double gap_rolling_ref[3] = {0.06858666958093279, 0.05062425428544091, -0.05782673704478762};
    double rollingcontactmoment_ref[3] = {
        -0.1119618097864935, 0.09699534835316952, -0.02050722172996138};

    BACI_EXPECT_ITERABLE_NEAR(gap_rolling, gap_rolling_ref, 3, 1.0e-12);

    for (int i = 0; i < 3; ++i)
      EXPECT_NEAR(rollingcontactmoment[i], rollingcontactmoment_ref[i], 1.0e-12);

    EXPECT_FALSE(stick_rolling);
  }

  TEST_F(DEMContactRollingCoulombTest, RollingContactMomentSlip)
  {
    double gap_rolling[3] = {0.0};
    gap_rolling[0] = 0.1;
    gap_rolling[1] = 0.05;
    gap_rolling[2] = -0.25;

    bool stick_rolling = false;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double v_rel_rolling[3] = {0.0};
    v_rel_rolling[0] = -0.03;
    v_rel_rolling[1] = 0.1;
    v_rel_rolling[2] = 0.12;

    const double m_eff = 2.5;
    const double r_eff = 0.5;
    const double normalcontactforce = 1.5;

    double rollingcontactmoment[3] = {0.0};
    contactrolling_->RollingContactMoment(gap_rolling, stick_rolling, e_ji, v_rel_rolling, m_eff,
        r_eff, mu_rolling_, normalcontactforce, rollingcontactmoment);

    double gap_rolling_ref[3] = {0.06858666958093279, 0.05062425428544091, -0.05782673704478762};
    double rollingcontactmoment_ref[3] = {
        -0.1119618097864935, 0.09699534835316952, -0.02050722172996138};

    BACI_EXPECT_ITERABLE_NEAR(gap_rolling, gap_rolling_ref, 3, 1.0e-12);

    for (int i = 0; i < 3; ++i)
      EXPECT_NEAR(rollingcontactmoment[i], rollingcontactmoment_ref[i], 1.0e-12);

    EXPECT_FALSE(stick_rolling);
  }

  TEST_F(DEMContactRollingCoulombTest, RollingPotentialEnergy)
  {
    double gap_rolling[3] = {0.0};
    gap_rolling[0] = 0.1;
    gap_rolling[1] = 0.05;
    gap_rolling[2] = -0.25;

    const double k_rolling = (1.0 - nue_) / (1.0 - 0.5 * nue_) * k_normal_per_radius_;
    const double rollingpotentialenergy_ref =
        0.5 * k_rolling * PARTICLEINTERACTION::UTILS::VecDot(gap_rolling, gap_rolling);

    double rollingpotentialenergy = 0.0;
    contactrolling_->RollingPotentialEnergy(gap_rolling, rollingpotentialenergy);

    EXPECT_NEAR(rollingpotentialenergy, rollingpotentialenergy_ref, 1.0e-12);
  }
}  // namespace
