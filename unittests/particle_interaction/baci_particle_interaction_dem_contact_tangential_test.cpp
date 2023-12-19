/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for tangential contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_particle_interaction_dem_contact_tangential.H"

#include "baci_inpar_validparameters.H"
#include "baci_particle_interaction_utils.H"
#include "baci_unittest_utils_assertions_test.H"

namespace
{
  using namespace BACI;

  class DEMContactTangentialLinearSpringDampTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp> contacttangential_;

    const double e_ = 0.8;
    const double nue_ = 0.4;
    const double mu_tangential_ = 0.2;

    const double k_normal_per_radius_ = 4.0;

    DEMContactTangentialLinearSpringDampTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("COEFF_RESTITUTION", e_);
      params_dem.set("POISSON_RATIO", nue_);
      params_dem.set("FRICT_COEFF_TANG", mu_tangential_);

      // create tangential contact handler
      contacttangential_ =
          std::make_unique<PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp>(params_dem);

      // init tangential contact handler
      contacttangential_->Init();

      // setup tangential contact handler
      contacttangential_->Setup(k_normal_per_radius_);
    }
    // note: the public functions Init() and Setup() of class DEMContactTangentialLinearSpringDamp
    // are called in the constructor and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactTangentialLinearSpringDampTest, TangentialContactForceStick)
  {
    double gap_tangential[3] = {0.0};
    gap_tangential[0] = 0.1;
    gap_tangential[1] = 0.05;
    gap_tangential[2] = -0.25;

    bool stick_tangential = true;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double vel_rel_tangential[3] = {0.0};
    vel_rel_tangential[0] = -0.03;
    vel_rel_tangential[1] = 0.1;
    vel_rel_tangential[2] = 0.12;

    const double m_eff = 2.5;
    const double normalcontactforce = 2.5e2;

    double tangentialcontactforce[3];
    contacttangential_->TangentialContactForce(gap_tangential, stick_tangential, e_ji,
        vel_rel_tangential, m_eff, mu_tangential_, normalcontactforce, tangentialcontactforce);

    double gap_tangential_ref[3] = {0.17923102017884, 0.16378007016342, -0.12669779012642};
    double tangentialcontactforce_ref[3] = {-0.52425016124507, -0.53614987479507, 0.32632177321349};

    BACI_EXPECT_ITERABLE_NEAR(gap_tangential, gap_tangential_ref, 3, 1.0e-14);

    for (int i = 0; i < 3; ++i)
      EXPECT_NEAR(tangentialcontactforce[i], tangentialcontactforce_ref[i], 1.0e-14);

    EXPECT_TRUE(stick_tangential);
  }

  TEST_F(DEMContactTangentialLinearSpringDampTest, TangentialContactForceSlip)
  {
    double gap_tangential[3] = {0.0};
    gap_tangential[0] = 0.1;
    gap_tangential[1] = 0.05;
    gap_tangential[2] = -0.25;

    bool stick_tangential = false;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double vel_rel_tangential[3] = {0.0};
    vel_rel_tangential[0] = -0.03;
    vel_rel_tangential[1] = 0.1;
    vel_rel_tangential[2] = 0.12;

    const double m_eff = 2.5;
    const double normalcontactforce = 1.5;

    double tangentialcontactforce[3];
    contacttangential_->TangentialContactForce(gap_tangential, stick_tangential, e_ji,
        vel_rel_tangential, m_eff, mu_tangential_, normalcontactforce, tangentialcontactforce);

    double gap_tangential_ref[3] = {0.068586669580933, 0.050624254285441, -0.057826737044788};
    double tangentialcontactforce_ref[3] = {-0.19231710945136, -0.19668242716113, 0.11970861396859};

    BACI_EXPECT_ITERABLE_NEAR(gap_tangential, gap_tangential_ref, 3, 1.0e-14);

    for (int i = 0; i < 3; ++i)
      EXPECT_NEAR(tangentialcontactforce[i], tangentialcontactforce_ref[i], 1.0e-14);

    EXPECT_FALSE(stick_tangential);
  }

  TEST_F(DEMContactTangentialLinearSpringDampTest, TangentialPotentialEnergy)
  {
    double gap_tangential[3] = {0.0};
    gap_tangential[0] = 0.1;
    gap_tangential[1] = 0.05;
    gap_tangential[2] = -0.25;

    const double k_tangential = (1.0 - nue_) / (1.0 - 0.5 * nue_) * k_normal_per_radius_;
    const double tangentialpotentialenergy_ref =
        0.5 * k_tangential * PARTICLEINTERACTION::UTILS::VecDot(gap_tangential, gap_tangential);

    double tangentialpotentialenergy = 0.0;
    contacttangential_->TangentialPotentialEnergy(gap_tangential, tangentialpotentialenergy);

    EXPECT_NEAR(tangentialpotentialenergy, tangentialpotentialenergy_ref, 1.0e-12);
  }
}  // namespace
