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
 * @file 	continuum_integration_1st_ck.h
 * @brief 	Here, we define the algorithm classes for continuum dynamics within the body.
 * @details CK and SYCL version.
 * @author	Shuang Li,Xiangyu Hu and Shuaihao Zhang
 */
#ifndef CONTINUUM_SEGREGATION_PARAMETERS_H
#define CONTINUUM_SEGREGATION_PARAMETERS_H

#include "constraint_dynamics.h"
#include "acoustic_step_1st_half.h"
#include "general_continuum.h"
#include "general_continuum.hpp"
namespace SPH
{
namespace continuum_dynamics
{

template <class BaseInteractionType>
class SegregationParametersUpdateBase : public BaseInteractionType
{

  public:
    template <class DynamicsIdentifier>
    explicit SegregationParametersUpdateBase(DynamicsIdentifier &identifier);
    virtual ~SegregationParametersUpdateBase(){};

  protected:
    BaseMaterial &material_;
    DiscreteVariable<Real> *dv_Vol_, *dv_seg_phi_;
    DiscreteVariable<Vecd> *dv_seg_n_;
    DiscreteVariable<int> *dv_indicator_;
};

template <typename...>
class UpdateSegNormalDirection;

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class UpdateSegNormalDirection<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public SegregationParametersUpdateBase<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = SegregationParametersUpdateBase<Interaction<Inner<Parameters...>>>;

  public:
    explicit UpdateSegNormalDirection(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~UpdateSegNormalDirection(){};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *seg_n_;
        Real *seg_phi_;
    };

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        KernelCorrectionType correction_;
        Real *Vol_;
        Vecd *seg_n_;
        Real *seg_phi_;
        int *indicator_;

        Real W0_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *seg_n_;
        Real *seg_phi_;
        int *indicator_;
    };

  protected:
    KernelCorrectionType correction_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class UpdateSegNormalDirection<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public SegregationParametersUpdateBase<Interaction<Contact<Wall, Parameters...>>>
{
    using BaseInteraction = SegregationParametersUpdateBase<Interaction<Contact<Wall, Parameters...>>>;

  public:
    explicit UpdateSegNormalDirection(Relation<Contact<Parameters...>> &wall_contact_relation);
    virtual ~UpdateSegNormalDirection(){};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        KernelCorrectionType correction_;
        
        Real *wall_Vol_;
        Vecd *seg_n_;
        Real *seg_phi_;
        int *indicator_;
    };
  protected:
    KernelCorrectionType correction_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class UpdateSegNormalDirection<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public SegregationParametersUpdateBase<Interaction<Contact<Soil, Parameters...>>>
{
    using BaseInteraction = SegregationParametersUpdateBase<Interaction<Contact<Soil, Parameters...>>>;

  public:
    explicit UpdateSegNormalDirection(Relation<Contact<Parameters...>> &wall_contact_relation);
    virtual ~UpdateSegNormalDirection(){};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        KernelCorrectionType correction_;
        
        Real *soil_Vol_;
        Vecd *seg_n_;
    };
  protected:
    KernelCorrectionType correction_;
};


using UpdateSegNormalDirectionRiemann =
    UpdateSegNormalDirection<Inner<OneLevel, AcousticRiemannSolver, NoKernelCorrection>,
                        Contact<Wall, AcousticRiemannSolver, NoKernelCorrection>>;
using WallUpdateSegNormalDirectionRiemann =
    UpdateSegNormalDirection<Inner<OneLevel, AcousticRiemannSolver, NoKernelCorrection>,
                        Contact<Soil, AcousticRiemannSolver, NoKernelCorrection>>;
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_SEGREGATION_PARAMETERS_HPP