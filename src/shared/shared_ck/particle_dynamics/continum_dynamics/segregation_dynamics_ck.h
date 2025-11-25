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
 * @file 	stress_diffusion.h
 * @brief 	Here, we define the ck_version for stress diffusion. 
 * @details Refer to Zhu et al(2023).
 * @author	Shuang Li, Xiangyu Hu
 */

#ifndef SEGREGATION_DYNAMICS_CK_H
#define SEGREGATION_DYNAMICS_CK_H

#include "base_continuum_dynamics.h"
#include "constraint_dynamics.h"
#include "fluid_integration.hpp"
#include "general_continuum.h"
#include "general_continuum.hpp"
#include "continuum_integration_1st_ck.h"
#include "continuum_integration_1st_ck.hpp"
#include "base_general_dynamics.h"
namespace SPH
{
namespace continuum_dynamics
{

class BaseSegregationParametersCK
{
  public:
    BaseSegregationParametersCK(BaseParticles *particles){};
    virtual ~BaseSegregationParametersCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0)
        {
        };

      protected:

    };

  protected:

};




class SegregationParametersCK : public LocalDynamics, public BaseSegregationParametersCK
{
  public:
    SegregationParametersCK(SPHBody &sph_body);
    virtual ~SegregationParametersCK() {};

    class UpdateKernel : public BaseSegregationParametersCK::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0)
        {
            diffusivity_[index_i] = 1.01e-3;
        };

      protected:
        Real *diffusivity_, *segregation_rate_;
    };

  protected:
      DiscreteVariable<Real> *dv_diffusivity_, *dv_segregation_rate_, *dv_segregation_test_;
};


} // namespace continuum_dynamics
} // namespace SPH
#endif // SEGREGATION_DYNAMICS_CK_H