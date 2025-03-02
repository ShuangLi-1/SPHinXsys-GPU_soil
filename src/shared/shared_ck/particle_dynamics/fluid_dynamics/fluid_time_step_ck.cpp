#include "fluid_time_step_ck.h"

#include "viscosity.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
AcousticTimeStepCK::AcousticTimeStepCK(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
AcousticTimeStepCK::FinishDynamics::FinishDynamics(AcousticTimeStepCK &encloser)
    : h_min_(encloser.h_min_), acousticCFL_(encloser.acousticCFL_) {}
//=================================================================================================//
Real AcousticTimeStepCK::FinishDynamics::Result(Real reduced_value)
{
    // since the particle does not change its configuration in the acoustic time steps
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
AdvectionTimeStepCK::
    AdvectionTimeStepCK(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")) {}
//=================================================================================================//
AdvectionTimeStepCK::FinishDynamics::FinishDynamics(AdvectionTimeStepCK &encloser)
    : h_min_(encloser.h_min_), speed_ref_(encloser.speed_ref_),
      advectionCFL_(encloser.advectionCFL_) {}
//=================================================================================================//
Real AdvectionTimeStepCK::FinishDynamics::Result(Real reduced_value)
{
    return advectionCFL_ * h_min_ / (SMAX(std::sqrt(reduced_value), speed_ref_) + TinyReal);
}
//=================================================================================================//
AdvectionViscousTimeStepCK::AdvectionViscousTimeStepCK(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : AdvectionTimeStepCK(sph_body, U_ref, advectionCFL)
{
    BaseMaterial &material = particles_->getBaseMaterial();
    Viscosity &viscosity = DynamicCast<Viscosity>(this, material);
    Real viscous_speed = viscosity.ReferenceViscosity() / material.ReferenceDensity() / h_min_;
    speed_ref_ = SMAX(viscous_speed, speed_ref_);
}
//=================================================================================================//
AdvectionStepSetup::AdvectionStepSetup(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      dv_Vol_(particles_->getVariableByName<Real>("VolumetricMeasure")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_dpos_(particles_->registerStateVariableOnly<Vecd>("Displacement")) {}
//=================================================================================================//
AdvectionStepClose::AdvectionStepClose(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_dpos_(particles_->getVariableByName<Vecd>("Displacement")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
