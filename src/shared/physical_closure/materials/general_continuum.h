/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	general_continuum.h
 * @brief 	Describe the linear elastic, J2 plasticity, and Drucker-Prager's plastic model
 * @author	Shuaihao Zhang and Xiangyu Hu
 */

#ifndef GENERAL_CONTINUUM_H
#define GENERAL_CONTINUUM_H

#include "weakly_compressible_fluid.h"

namespace SPH
{
class GeneralContinuum : public WeaklyCompressibleFluid
{
  protected:
    Real E_;                 /* Youngs or tensile modules  */
    Real G_;                 /* shear modules  */
    Real K_;                 /* bulk modules  */
    Real nu_;                /* Poisson ratio  */
    Real contact_stiffness_; /* contact-force stiffness related to bulk modulus*/
  public:
    explicit GeneralContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio)
        : WeaklyCompressibleFluid(rho0, c0), E_(0.0), G_(0.0), K_(0.0), nu_(0.0), contact_stiffness_(rho0_ * c0 * c0)
    {
        material_type_name_ = "GeneralContinuum";
        E_ = youngs_modulus;
        nu_ = poisson_ratio;
        G_ = getShearModulus(youngs_modulus, poisson_ratio);
        K_ = getBulkModulus(youngs_modulus, poisson_ratio);
        lambda0_ = getLambda(youngs_modulus, poisson_ratio);
    };
    virtual ~GeneralContinuum(){};

    Real lambda0_; /* first Lame parameter */
    Real getYoungsModulus() { return E_; };
    Real getPoissonRatio() { return nu_; };
    Real getDensity() { return rho0_; };
    Real getBulkModulus(Real youngs_modulus, Real poisson_ratio);
    Real getShearModulus(Real youngs_modulus, Real poisson_ratio);
    Real getLambda(Real youngs_modulus, Real poisson_ratio);

    Real ContactStiffness() { return contact_stiffness_; };

    virtual Matd ConstitutiveRelationShearStress(Matd &velocity_gradient, Matd &shear_stress);


    class GeneralContinuumKernel : public WeaklyCompressibleFluid::EosKernel
    {
      public:
        GeneralContinuumKernel(GeneralContinuum &encloser) : WeaklyCompressibleFluid::EosKernel(encloser),
        E_(encloser.E_), G_(encloser.G_),K_(encloser.K_),
        nu_(encloser.nu_),contact_stiffness_(encloser.contact_stiffness_),
        rho0_(encloser.rho0_){};

        inline Real getYoungsModulus() { return E_; };
        inline Real getPoissonRatio() { return nu_; };
        inline Real getDensity() { return rho0_; };
        inline Real getBulkModulus(Real youngs_modulus, Real poisson_ratio);
        inline Real getShearModulus(Real youngs_modulus, Real poisson_ratio);
        inline Real getLambda(Real youngs_modulus, Real poisson_ratio);
      protected:
        Real E_;                 /* Youngs or tensile modules  */
        Real G_;                 /* shear modules  */
        Real K_;                 /* bulk modules  */
        Real nu_;                /* Poisson ratio  */
        Real contact_stiffness_; /* contact-force stiffness related to bulk modulus*/
        Real rho0_; /* contact-force stiffness related to bulk modulus*/
    };
};

class PlasticContinuum : public GeneralContinuum
{
  protected:
    Real c_;                            /* cohesion  */
    Real phi_;                          /* friction angle  */
    Real psi_;                          /* dilatancy angle  */
    Real alpha_phi_;                    /* Drucker-Prager's constants */
    Real k_c_;                          /* Drucker-Prager's constants */
    Real miu_s_, miu_d_, I0_;                 /* Miu(I) constants */
    Real d_min_, d_max_;
    const Real stress_dimension_ = 3.0; /* plain strain condition */
  public:
    explicit PlasticContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio, Real friction_angle,
     Real cohesion = 0, Real dilatancy = 0, Real viscosity_s = 0.384, Real viscosity_d = 0.65, Real d_min = 0.002, Real d_max = 0.002)
        : GeneralContinuum(rho0, c0, youngs_modulus, poisson_ratio),
          c_(cohesion), phi_(friction_angle), psi_(dilatancy), alpha_phi_(0.0), k_c_(0.0),
          miu_s_(viscosity_s), miu_d_(viscosity_d), d_min_(d_min), d_max_(d_max)
    {
        I0_ = 0.279;
        material_type_name_ = "PlasticContinuum";
        alpha_phi_ = getDPConstantsA(friction_angle);
        k_c_ = getDPConstantsK(cohesion, friction_angle);
    };
    /*7 Parameters*/
    explicit PlasticContinuum(ConstructArgs<Real, Real, Real, Real, Real, Real, Real> args)
    : PlasticContinuum(std::get<0>(args), std::get<1>(args), std::get<2>(args), std::get<3>(args), 
    std::get<4>(args), std::get<5>(args), std::get<6>(args)) {};
    /*9 Parameters*/
    explicit PlasticContinuum(ConstructArgs<Real, Real, Real, Real, Real, Real, Real, Real, Real> args)
    : PlasticContinuum(std::get<0>(args), std::get<1>(args), std::get<2>(args), std::get<3>(args), 
    std::get<4>(args), std::get<5>(args), std::get<6>(args), std::get<7>(args), std::get<8>(args)) {};
    /*11 Parameters*/
    explicit PlasticContinuum(ConstructArgs<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real> args)
    : PlasticContinuum(std::get<0>(args), std::get<1>(args), std::get<2>(args), std::get<3>(args), 
    std::get<4>(args), std::get<5>(args), std::get<6>(args), std::get<7>(args), std::get<8>(args),
    std::get<9>(args),std::get<10>(args)) {};
    virtual ~PlasticContinuum(){};

    Real getDPConstantsA(Real friction_angle);
    Real getDPConstantsK(Real cohesion, Real friction_angle);
    Real getDPConstantsA_WithMiu(Real miu);
    Real getDPConstantsK_WithMiu(Real miu);
    Real getDPConstantsK_WithMiu(Real cohesion, Real miu);
    Real getFrictionAngle() { return phi_; };
    Real getCohesion() { return c_; };
    Real getStaticViscosity() {return miu_s_; };
    Real getDynamicViscosity() {return miu_d_; };
    Real getMaterialConstant() {return I0_; };
    Real getMinDiameter() {return d_min_; };
    Real getMaxDiameter() {return d_max_; };

    virtual Mat3d ConstitutiveRelation(Mat3d &velocity_gradient, Mat3d &stress_tensor);
    virtual Mat3d ReturnMapping(Mat3d &stress_tensor);
    virtual Mat3d ConstitutiveRelation_withMiuI(Mat3d &velocity_gradient, Mat3d &stress_tensor, Real alpha_phi_i, Real k_c_i);
    virtual Mat3d ReturnMapping_withMiuI(Mat3d &stress_tensor, Real alpha_phi_i, Real k_c_i);
    virtual Mat3d NEW_ReturnMapping_withMiuI(Mat3d &stress_tensor, Real alpha_phi_i, Real k_c_i);
    virtual Real getViscosity(Real p, Real equivalentShearStrainRate, Real yita_0 = 0.0)
    {
      if(equivalentShearStrainRate > TinyReal)
        return yita_0 + (c_ + p* tan(phi_) )/equivalentShearStrainRate;
      else
        return 0.0;
    }
    inline Real getViscosityMiuI(Real inertial_number_i, Real miu_s=0.365, Real miu_d=0.572)
    {
        Real inertial_number_0=0.279;
        Real miu = miu_s + (miu_d - miu_s)/(inertial_number_i/inertial_number_0 + 1.0);
        return miu;
    }



    class PlasticKernel: public GeneralContinuum::GeneralContinuumKernel
    {
      public:

        PlasticKernel(PlasticContinuum &encloser) : GeneralContinuum::GeneralContinuumKernel(encloser),
        c_(encloser.c_),phi_(encloser.phi_),
        psi_(encloser.psi_),alpha_phi_(encloser.alpha_phi_),k_c_(encloser.k_c_){};

        inline Real getDPConstantsA(Real friction_angle);
        inline Mat3d ConstitutiveRelation(Mat3d &velocity_gradient, Mat3d &stress_tensor);  
        inline Mat3d ReturnMapping(Mat3d &stress_tensor);
        inline Real getFrictionAngle() { return phi_; };
        inline Real getViscosity(Real p, Real equivalentShearStrainRate)
        {
          Real yita_0 = 0.0;
          Real cohesion = 0.0;
          Real tan_fai = 0.3;
          if(equivalentShearStrainRate > TinyReal)
            return yita_0 + (cohesion + p*tan_fai)/equivalentShearStrainRate;
          else
            return 0.0;
        }

      protected:
          Real c_;                            /* cohesion  */
          Real phi_;                          /* friction angle  */
          Real psi_;                          /* dilatancy angle  */
          Real alpha_phi_;                    /* Drucker-Prager's constants */
          Real k_c_;                          /* Drucker-Prager's constants */
          Real stress_dimension_ = 3.0; /* plain strain condition */   //Temporarily cancel const --need to check


    };
};

class J2Plasticity : public GeneralContinuum
{
  protected:
    Real yield_stress_;
    Real hardening_modulus_;
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);

  public:
    explicit J2Plasticity(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio, Real yield_stress, Real hardening_modulus = 0.0)
        : GeneralContinuum(rho0, c0, youngs_modulus, poisson_ratio),
          yield_stress_(yield_stress), hardening_modulus_(hardening_modulus)
    {
        material_type_name_ = "J2Plasticity";
    };
    virtual ~J2Plasticity(){};

    Real YieldStress() { return yield_stress_; };
    Real HardeningModulus() { return hardening_modulus_; };

    virtual Matd ConstitutiveRelationShearStress(Matd &velocity_gradient, Matd &shear_stress, Real &hardening_factor);
    virtual Matd ReturnMappingShearStress(Matd &shear_stress,  Real &hardening_factor);
    virtual Real ScalePenaltyForce(Matd &shear_stress, Real &hardening_factor);
    virtual Real HardeningFactorRate(const Matd &shear_stress, Real &hardening_factor);
};
} // namespace SPH
#endif // GENERAL_CONTINUUM_H