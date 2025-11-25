#pragma once

#include "segregation_shape.h"

namespace SPH
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
SegregationNormalAndDegeneration<Base, DataDelegationType>::SegregationNormalAndDegeneration(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      seg_n_(this->particles_->template registerStateVariable<Vecd>("SegNormalDirection")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      seg_phi_(this->particles_->template registerStateVariable<Real>("SegPhi")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      gama_(this->particles_->template registerStateVariable<Real>("ShearRate")),
      indicator_(this->particles_->template getVariableDataByName<int>("Indicator")),
      velocity_gradient_(this->particles_->template registerStateVariable<Matd>("VelocityGradient")),
      W0_(this->sph_body_.getSPHAdaptation().getKernel()->W0(ZeroVecd)) {}
//=================================================================================================//


} // namespace SPH
