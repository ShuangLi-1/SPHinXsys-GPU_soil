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
#ifndef NONLOCAL_H
#define NONLOCAL_H

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
class NonLocalPlasticIntegration1stHalf;

template <class RiemannSolverType>
class NonLocalPlasticIntegration1stHalf<Inner<>, RiemannSolverType>
    : public BasePlasticIntegration<DataDelegateInner>
{
  public:
    explicit NonLocalPlasticIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~NonLocalPlasticIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
    virtual Vecd computeNonConservativeForce(size_t index_i);

  protected:
    RiemannSolverType riemann_solver_;
};
using NonLocalPlasticIntegration1stHalfInnerNoRiemann = NonLocalPlasticIntegration1stHalf<Inner<>, NoRiemannSolver>;
using NonLocalPlasticIntegration1stHalfInnerRiemann = NonLocalPlasticIntegration1stHalf<Inner<>, AcousticRiemannSolver>;


template <class RiemannSolverType>
class NonLocalPlasticIntegration1stHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit NonLocalPlasticIntegration1stHalf(BaseContactRelation &wall_contact_relation);
    virtual ~NonLocalPlasticIntegration1stHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
    virtual Vecd computeNonConservativeForce(size_t index_i);

  protected:
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType>
using NonLocalPlasticIntegration1stHalfWithWall = ComplexInteraction<NonLocalPlasticIntegration1stHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;
using NonLocalPlasticIntegration1stHalfWithWallNoRiemann = NonLocalPlasticIntegration1stHalfWithWall<NoRiemannSolver>;
using NonLocalPlasticIntegration1stHalfWithWallRiemann = NonLocalPlasticIntegration1stHalfWithWall<AcousticRiemannSolver>;

template <typename... InteractionTypes>
class NonLocalPlasticIntegration2ndHalf;

template <class RiemannSolverType>
class NonLocalPlasticIntegration2ndHalf<Inner<>, RiemannSolverType>
    : public BasePlasticIntegration<DataDelegateInner>
{
  public:
    explicit NonLocalPlasticIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~NonLocalPlasticIntegration2ndHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    Real *Vol_, *mass_;
    Real *concentration_;
};
using NonLocalPlasticIntegration2ndHalfInnerNoRiemann = NonLocalPlasticIntegration2ndHalf<Inner<>, NoRiemannSolver>;
using NonLocalPlasticIntegration2ndHalfInnerRiemann = NonLocalPlasticIntegration2ndHalf<Inner<>, AcousticRiemannSolver>;

template <class RiemannSolverType>
class NonLocalPlasticIntegration2ndHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit NonLocalPlasticIntegration2ndHalf(BaseContactRelation &wall_contact_relation);
    virtual ~NonLocalPlasticIntegration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    Real *concentration_;
};

template <class RiemannSolverType>
using NonLocalPlasticIntegration2ndHalfWithWall = ComplexInteraction<NonLocalPlasticIntegration2ndHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;
using NonLocalPlasticIntegration2ndHalfWithWallNoRiemann = NonLocalPlasticIntegration2ndHalfWithWall<NoRiemannSolver>;
using NonLocalPlasticIntegration2ndHalfWithWallRiemann = NonLocalPlasticIntegration2ndHalfWithWall<AcousticRiemannSolver>;


/**
 * @class UpdateViscosityAndParticleTemp
 */
class UpdateNonLocalParameters : public LocalDynamics
{
  public:
    explicit UpdateNonLocalParameters(SPHBody &sph_body, Real scale_factor = 1.0);
    virtual ~UpdateNonLocalParameters(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    PlasticContinuum &plastic_continuum_;
    Real *test_;
    Real *concentration_, *gama_, *segregation_rate_, *diffusion_rate_, *inertial_num_, *p_, *rho_;
    Real *alpha_each_, *Kc_each_;
    Real *friction_, *fluidity_, *local_fluidity_rate_, *nonlocal_fluidity_rate_;
    Real scale_factor_, d_min_, d_max_, x_d_, x_v_;
    Real d_mean_, A_, beta_, R_, Fai_, C_, epsilon_, actual_height_, gravity_;
    Real rho0_, cohesion_, miu_s_, miu_d_, I0_, t0_;
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // NONLOCAL_H