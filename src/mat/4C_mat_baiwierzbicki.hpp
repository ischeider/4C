// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_BAIWIERZBICKI_HPP
#define FOUR_C_MAT_BAIWIERZBICKI_HPP

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core
{
  namespace Utils
  {
    class FunctionOfAnything;
  }
}  // namespace Core

namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for the Modified Bai-Wierzbicki (MBW) damage model
    //!
    //! Large-strain (logarithmic / Hencky strain, additive split assuming small elastic
    //! strains, nonlinear kinematics after Miehe) von Mises plasticity coupled with a
    //! hybrid continuum-damage-mechanics ductile damage model according to
    //! Lian et al., Int. J. Damage Mech. 24(3), 2015.
    class BaiWierzbicki : public Core::Mat::PAR::Parameter
    {
     public:
      //! standard constructor
      BaiWierzbicki(const Core::Mat::PAR::Parameter::Data& matdata);

      //! @name material parameters
      //@{

      //! mass density
      const double density_;
      //! Young's modulus
      const double youngs_;
      //! Poisson's ratio
      const double poissonratio_;
      //! initial (constant) yield stress, used if HARDENING_FUNC == 0
      const double yield_;
      //! function ID for the isotropic hardening curve sigma_y(epsp)
      const int functionID_hardening_;
      //! function ID for the damage initiation locus epsp_i(eta, theta)
      const int functionID_dil_;
      //! specific fracture energy G_f
      const double failenergy_;
      //! critical damage value D_cr at which the load carrying capacity is lost
      const double dcrit_;
      //! integrity limit omega = 1 - D below which the Gauss point is considered failed
      const double omega_failure_;
      //! viscosity for the Perzyna viscous regularization
      const double visc_;
      //! rate dependency exponent for the Perzyna viscous regularization
      const double rate_dependency_;
      //! characteristic length; if 0, the characteristic Gauss point size is used
      const double charlength_;
      //! maximum number of local Newton-Raphson iterations
      const int max_iterations_;
      //! tolerance for the local Newton-Raphson iteration
      const double tolerance_nr_;

      //@}

      //! create material instance of matching type with my parameters
      std::shared_ptr<Core::Mat::Material> create_material() override;

    };  // class BaiWierzbicki

  }  // namespace PAR


  class BaiWierzbickiType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "BaiWierzbickiType"; }

    static BaiWierzbickiType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static BaiWierzbickiType instance_;

  };  // class BaiWierzbickiType

  /*----------------------------------------------------------------------*/
  //! wrapper for the Modified Bai-Wierzbicki ductile damage material
  class BaiWierzbicki : public So3Material
  {
   public:
    //! construct empty material object
    BaiWierzbicki();

    //! construct the material object given material parameters
    explicit BaiWierzbicki(Mat::PAR::BaiWierzbicki* params);

    //! @name Packing and Unpacking

    int unique_par_object_id() const override
    {
      return BaiWierzbickiType::instance().unique_par_object_id();
    }

    void pack(Core::Communication::PackBuffer& data) const override;

    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    //@}

    //! @name Access methods

    //! material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_baiwierzbicki;
    }

    /// check if element kinematics and material kinematics are compatible
    void valid_kinematics(Solid::KinemType kinem) override
    {
      if (!(kinem == Solid::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    //! return copy of this material object
    std::shared_ptr<Core::Mat::Material> clone() const override
    {
      return std::make_shared<BaiWierzbicki>(*this);
    }

    //! density
    double density() const override { return params_->density_; }

    //! return quick accessible material parameter data
    Core::Mat::PAR::Parameter* parameter() const override { return params_; }

    //! check if history variables are already initialized
    bool initialized() const { return (isinit_ and !accplstraincurr_.empty()); }

    //! return accumulated (equivalent) plastic strain at a Gauss point (current step)
    double accumulated_strain(int gp) const { return accplstraincurr_.at(gp); }

    //! return the damage variable at a Gauss point (current step)
    double damage(int gp) const { return damagecurr_.at(gp); }

    //! return the accumulated plastic dissipated energy at a Gauss point (current step)
    double plastic_dissipated_energy(int gp) const { return plastenergycurr_.at(gp); }

    //! return the accumulated artificial (viscous regularization) energy at a Gauss point
    double viscous_artificial_energy(int gp) const { return viscenergycurr_.at(gp); }

    //! return names of visualization data
    void vis_names(std::map<std::string, int>& names) const override;

    //! return visualization data
    bool vis_data(
        const std::string& name, std::vector<double>& data, int numgp, int eleID) const override;

    //! return names of visualization data available for direct VTK output
    void register_output_data_names(
        std::unordered_map<std::string, int>& names_and_size) const override;

    //! return visualization data for direct VTK output
    bool evaluate_output_data(
        const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const override;

    //@}

    //! @name Evaluation methods

    //! initialise internal stress variables
    void setup(int numgp, const Discret::Elements::Fibers& fibers,
        const std::optional<Discret::Elements::CoordinateSystem>& coord_system) override;

    //! update internal stress variables
    void update() override;

    //! evaluate material law
    void evaluate(const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
        const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
        Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
        Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID) override;

   private:
    //! my material parameters
    Mat::PAR::BaiWierzbicki* params_;

    //! inverse plastic right Cauchy-Green tensor C_p^{-1}
    std::vector<Core::LinAlg::SymmetricTensor<double, 3, 3>> invplrcglast_;
    std::vector<Core::LinAlg::SymmetricTensor<double, 3, 3>> invplrcgcurr_;

    //! accumulated (equivalent) plastic strain
    std::vector<double> accplstrainlast_;
    std::vector<double> accplstraincurr_;

    //! isotropic damage variable D
    std::vector<double> damagelast_;
    std::vector<double> damagecurr_;

    //! damage initiation indicator w = int deps_p / epsp_i (initiation at w >= 1)
    std::vector<double> daminitindiclast_;
    std::vector<double> daminitindiccurr_;

    //! yield stress frozen at damage initiation (sigma_y0); negative means not yet initiated
    std::vector<double> yieldinitlast_;
    std::vector<double> yieldinitcurr_;

    //! accumulated plastic dissipated energy density
    std::vector<double> plastenergylast_;
    std::vector<double> plastenergycurr_;

    //! accumulated artificial (viscous regularization) energy density
    std::vector<double> viscenergylast_;
    std::vector<double> viscenergycurr_;

    //! Gauss point failure flag (0 = intact, 1 = failed)
    std::vector<double> failedlast_;
    std::vector<double> failedcurr_;

    //! characteristic length associated with each Gauss point (constant)
    std::vector<double> charlength_;

    //! output-only quantities at t_{n+1} (recomputed at each evaluation)
    std::vector<double> triaxiality_;
    std::vector<double> lodeparam_;
    std::vector<double> dilvalue_;

    const Core::Utils::FunctionOfAnything* hardening_function_{nullptr};
    const Core::Utils::FunctionOfAnything* dil_function_{nullptr};

    //! indicator if setup routine has been called
    bool isinit_{false};
  };  // class BaiWierzbicki
}  // namespace Mat

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
