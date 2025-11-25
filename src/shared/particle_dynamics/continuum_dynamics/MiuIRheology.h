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
 * @file 	continuum_integration.h
 * @brief 	Here, we define the algorithm classes for continuum dynamics within the body.
 * @details We consider here weakly compressible assumption to model elastic and
 * 			plastic materials with the updated Lagrangian framework.
 * @author	Shuaihao Zhang and Xiangyu Hu
 */
#ifndef MIUIRHELOGY_H
#define MIUIRHELOGY_H

#include "base_continuum_dynamics.h"
#include "constraint_dynamics.h"
#include "fluid_integration.hpp"
#include "general_continuum.h"
#include "general_continuum.hpp"
#include "continuum_integration.h"
#include "continuum_integration.hpp"
namespace SPH
{
namespace continuum_dynamics
{

template <typename... InteractionTypes>
class MiuPlasticIntegration1stHalf;

template <class RiemannSolverType>
class MiuPlasticIntegration1stHalf<Inner<>, RiemannSolverType>
    : public BasePlasticIntegration<DataDelegateInner>
{
  public:
    explicit MiuPlasticIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~MiuPlasticIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
    virtual Vecd computeNonConservativeForce(size_t index_i);

  protected:
    RiemannSolverType riemann_solver_;
};
using MiuPlasticIntegration1stHalfInnerNoRiemann = MiuPlasticIntegration1stHalf<Inner<>, NoRiemannSolver>;
using MiuPlasticIntegration1stHalfInnerRiemann = MiuPlasticIntegration1stHalf<Inner<>, AcousticRiemannSolver>;


template <class RiemannSolverType>
class MiuPlasticIntegration1stHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit MiuPlasticIntegration1stHalf(BaseContactRelation &wall_contact_relation);
    virtual ~MiuPlasticIntegration1stHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
    virtual Vecd computeNonConservativeForce(size_t index_i);

  protected:
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType>
using MiuPlasticIntegration1stHalfWithWall = ComplexInteraction<MiuPlasticIntegration1stHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;
using MiuPlasticIntegration1stHalfWithWallNoRiemann = MiuPlasticIntegration1stHalfWithWall<NoRiemannSolver>;
using MiuPlasticIntegration1stHalfWithWallRiemann = MiuPlasticIntegration1stHalfWithWall<AcousticRiemannSolver>;

template <typename... InteractionTypes>
class MiuPlasticIntegration2ndHalf;

template <class RiemannSolverType>
class MiuPlasticIntegration2ndHalf<Inner<>, RiemannSolverType>
    : public BasePlasticIntegration<DataDelegateInner>
{
  public:
    explicit MiuPlasticIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~MiuPlasticIntegration2ndHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    Real *Vol_, *mass_;
};
using MiuPlasticIntegration2ndHalfInnerNoRiemann = MiuPlasticIntegration2ndHalf<Inner<>, NoRiemannSolver>;
using MiuPlasticIntegration2ndHalfInnerRiemann = MiuPlasticIntegration2ndHalf<Inner<>, AcousticRiemannSolver>;

template <class RiemannSolverType>
class MiuPlasticIntegration2ndHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit MiuPlasticIntegration2ndHalf(BaseContactRelation &wall_contact_relation);
    virtual ~MiuPlasticIntegration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType>
using MiuPlasticIntegration2ndHalfWithWall = ComplexInteraction<MiuPlasticIntegration2ndHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;
using MiuPlasticIntegration2ndHalfWithWallNoRiemann = MiuPlasticIntegration2ndHalfWithWall<NoRiemannSolver>;
using MiuPlasticIntegration2ndHalfWithWallRiemann = MiuPlasticIntegration2ndHalfWithWall<AcousticRiemannSolver>;

/**
 * @class UpdateViscosity
 */
class UpdateMiuIParameters : public LocalDynamics
{
  public:
    explicit UpdateMiuIParameters(SPHBody &sph_body, Real scale_factor = 1.0);
    virtual ~UpdateMiuIParameters(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    PlasticContinuum &plastic_continuum_;
    Real *concentration_, *gama_, *segregation_rate_, *diffusion_rate_, *inertial_num_, *p_, *rho_;
    Real *alpha_each_, *Kc_each_;
    Real scale_factor_, d_min_, d_max_, x_d_, x_v_;
    Real d_mean_, A_, beta_, R_, Fai_, C_, epsilon_, actual_height_, gravity_;
};

} // namespace continuum_dynamics
} // namespace SPH
#endif // MIUIRHELOGY_H