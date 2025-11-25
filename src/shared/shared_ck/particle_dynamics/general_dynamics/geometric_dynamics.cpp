#include "geometric_dynamics.h"

namespace SPH
{
//=================================================================================================//
NormalFromBodyShapeCK::NormalFromBodyShapeCK(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      initial_shape_(&sph_body.getInitialShape()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_n_(particles_->registerStateVariableOnly<Vecd>("NormalDirection")),
      dv_n0_(particles_->registerStateVariableOnly<Vecd>("InitialNormalDirection", dv_n_)),
      dv_phi_(particles_->registerStateVariableOnly<Real>("SignedDistance")),
      dv_phi0_(particles_->registerStateVariableOnly<Real>("InitialSignedDistance", dv_phi_)) {}
//=============================================================================================//
void NormalFromBodyShapeCK::UpdateKernel::update(size_t index_i, Real dt)
{
    Vecd normal_direction = initial_shape_->findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
    Real signed_distance = initial_shape_->findSignedDistance(pos_[index_i]);
    phi_[index_i] = signed_distance;
    phi0_[index_i] = signed_distance;
}
//=============================================================================================//
SurfaceIndicationFromBodyShape::SurfaceIndicationFromBodyShape(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      initial_shape_(&sph_body.getInitialShape()),
      spacing_ref_(sph_body.getSPHAdaptation().ReferenceSpacing()),
      dv_indicator_(particles_->registerStateVariableOnly<int>("SurfaceIndicator")),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")) {}
//=============================================================================================//
void SurfaceIndicationFromBodyShape::UpdateKernel::update(size_t index_i, Real dt)
{
    Real signed_distance = initial_shape_->findSignedDistance(pos_[index_i]);
    indicator_[index_i] = signed_distance > -spacing_ref_ ? 1 : 0;
}
//=================================================================================================//
NormalFromParticlesCK::NormalFromParticlesCK(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation), 
      initial_shape_(&sph_body_.getInitialShape()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_n_(particles_->registerStateVariableOnly<Vecd>("NormalDirection")),
      dv_n0_(particles_->registerStateVariableOnly<Vecd>("InitialNormalDirection", dv_n_)),
      dv_phi_(particles_->registerStateVariableOnly<Real>("SignedDistance")),
      dv_phi0_(particles_->registerStateVariableOnly<Real>("InitialSignedDistance", dv_phi_)),
      dv_Vol_(particles_->getVariableByName<Real>("Volume")) {}
      //=============================================================================================//
void NormalFromParticlesCK::InteractKernel::interact(size_t index_i, Real dt)
{
    // Vecd normal_direction = ZeroData<Vecd>::value;
    // const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    // for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    // {
    //     size_t index_j = inner_neighborhood.j_[n];
    //     normal_direction -= inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
    // }
    // normal_direction = normal_direction / (normal_direction.norm() + TinyReal);
    // n_[index_i] = normal_direction;
    // n0_[index_i] = normal_direction;
    // Real signed_distance = initial_shape_.findSignedDistance(pos_[index_i]);
    // phi_[index_i] = signed_distance;
    // phi0_[index_i] = signed_distance;
}
//=================================================================================================//
} // namespace SPH
