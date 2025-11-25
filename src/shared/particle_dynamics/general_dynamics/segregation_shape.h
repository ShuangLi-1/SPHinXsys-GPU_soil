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
 * @file    general_geometric.h
 * @brief   This is the particle dynamics applicable for all type bodies
 * @author	Xiangyu Hu
 */

#ifndef SEGREGATION_SHAPE_H
#define SEGREGATION_SHAPE_H

#include "base_general_dynamics.h"
#include "general_continuum.h"
#include "general_continuum.hpp"
namespace SPH
{
template <typename... InteractionTypes>
class SegregationNormalAndDegeneration;

template <class DataDelegationType>
class SegregationNormalAndDegeneration<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit SegregationNormalAndDegeneration(BaseRelationType &base_relation);
    virtual ~SegregationNormalAndDegeneration(){};

  protected:
    Vecd *seg_n_, *vel_;
    Real *seg_phi_, *Vol_, *gama_;
    int *indicator_;
    Matd *velocity_gradient_;
    Real W0_;
};


template <>
class SegregationNormalAndDegeneration<Inner<Base>> : public SegregationNormalAndDegeneration<Base, DataDelegateInner>
{
  public:
    explicit SegregationNormalAndDegeneration(BaseInnerRelation &inner_relation)
        : SegregationNormalAndDegeneration<Base, DataDelegateInner>(inner_relation){};
    virtual ~SegregationNormalAndDegeneration(){};
};


template <>
class SegregationNormalAndDegeneration<Inner<>> : public SegregationNormalAndDegeneration<Inner<Base>>
{
  public:
    explicit SegregationNormalAndDegeneration(BaseInnerRelation &inner_relation)
        : SegregationNormalAndDegeneration<Inner<Base>>(inner_relation){};
    virtual ~SegregationNormalAndDegeneration(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};


template <>
class SegregationNormalAndDegeneration<Contact<Base>> : public SegregationNormalAndDegeneration<Base, DataDelegateContact>
{
  public:
    explicit SegregationNormalAndDegeneration(BaseContactRelation &contact_relation);
    virtual ~SegregationNormalAndDegeneration(){};

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_;
};

template <>
class SegregationNormalAndDegeneration<Contact<Wall>> : public SegregationNormalAndDegeneration<Contact<Base>>
{
  public:
    explicit SegregationNormalAndDegeneration(BaseContactRelation &contact_relation)
        : SegregationNormalAndDegeneration<Contact<Base>>(contact_relation){};
    virtual ~SegregationNormalAndDegeneration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};


template <>
class SegregationNormalAndDegeneration<Contact<Soil>> : public SegregationNormalAndDegeneration<Contact<Base>>
{
  public:
    explicit SegregationNormalAndDegeneration(BaseContactRelation &contact_relation)
        : SegregationNormalAndDegeneration<Contact<Base>>(contact_relation){};
    virtual ~SegregationNormalAndDegeneration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseSegregationNormalandDegenerationComplex = ComplexInteraction<SegregationNormalAndDegeneration<InnerInteractionType, ContactInteractionTypes...>>;

using SoilSegAndPhiComplex = BaseSegregationNormalandDegenerationComplex<Inner<>, Contact<Wall>>;
using WallSegAndPhiComplex = BaseSegregationNormalandDegenerationComplex<Inner<>, Contact<Soil>>;

/**
 * @class NormalDirectionFromBodyShape
 * @brief normal direction at particles
 */
class UpdateDiffusivityAndSegregationRate : public LocalDynamics
{
  public:
    explicit UpdateDiffusivityAndSegregationRate(SPHBody &sph_body, Real scale_factor = 1.0);
    virtual ~UpdateDiffusivityAndSegregationRate(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    PlasticContinuum &plastic_continuum_;
    Real *test_;
    Real *concentration_, *gama_, *segregation_rate_, *diffusion_rate_, *inertial_num_, *p_, *rho_;
    Real *alpha_each_, *Kc_each_;
    Real scale_factor_, d_min_, d_max_, x_d_, x_v_;
    Real d_mean_, A_, beta_, R_, Fai_, C_, epsilon_, actual_height_, gravity_;
};

class UpdateDiffusivityAndSegregationRateGray : public LocalDynamics
{
  public:
    explicit UpdateDiffusivityAndSegregationRateGray(SPHBody &sph_body, Real height, Real scale_factor = 1.0);
    virtual ~UpdateDiffusivityAndSegregationRateGray(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *concentration_, *gama_, *segregation_rate_, *diffusion_rate_;
    Vecd *pos_;
    Real height_, period_, d_min_, d_max_, A_, beta_, R_, Fai_, C_, epsilon_, actual_height_;
};
class UpdateDiffusivityAndSegregationRateGray2023NonDim : public LocalDynamics
{
  public:
    explicit UpdateDiffusivityAndSegregationRateGray2023NonDim(SPHBody &sph_body, Real height, Real scale_factor = 1.0);
    virtual ~UpdateDiffusivityAndSegregationRateGray2023NonDim(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *concentration_, *gama_, *segregation_rate_, *diffusion_rate_;
    Vecd *pos_;
    Real height_, period_, d_min_, d_max_, d_mean_, A_, beta_, R_, Fai_, C_, epsilon_, actual_height_, U_;
};
class UpdateDiffusivityAndSegregationRateGray2021 : public LocalDynamics
{
  public:
    explicit UpdateDiffusivityAndSegregationRateGray2021(SPHBody &sph_body, Real height, Real scale_factor = 1.0);
    virtual ~UpdateDiffusivityAndSegregationRateGray2021(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *concentration_, *gama_, *segregation_rate_, *diffusion_rate_, *rho_, *p_;
    Vecd *pos_;
    Real height_, period_, d_min_, d_max_, d_mean_, A_, beta_, R_, Fai_, C_, epsilon_, actual_height_, gravity_;
};
} // namespace SPH
#endif // GENERAL_GEOMETRIC_H
