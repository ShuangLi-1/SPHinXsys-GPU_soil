#ifndef SEGREGATION_DYNAMICS_CK_HPP
#define SEGREGATION_DYNAMICS_CK_HPP

#include "segregation_dynamics_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
BaseSegregationParametersCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser){}
//=================================================================================================//
SegregationParametersCK::SegregationParametersCK(SPHBody &sph_body)
    : LocalDynamics(sph_body), BaseSegregationParametersCK(this->particles_),
    dv_diffusivity_(particles_->getVariableByName<Real>("SegregationDiffusivity")),
    dv_segregation_rate_(particles_->getVariableByName<Real>("SegregationRate")),
    dv_segregation_test_(particles_->getVariableByName<Real>("SegregationTest"))
    {
        particles_->addEvolvingVariable<Real>("SegregationDiffusivity");
        particles_->addEvolvingVariable<Real>("SegregationRate");
        particles_->addEvolvingVariable<Real>("SegregationTest");
    }
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
SegregationParametersCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseSegregationParametersCK::UpdateKernel(ex_policy, encloser),
    diffusivity_(encloser.dv_diffusivity_->DelegatedData(ex_policy)){}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH

#endif // SEGREGATION_DYNAMICS_CK_HPP