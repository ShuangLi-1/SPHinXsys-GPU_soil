#include "diffusion_segregation_reaction.h"
#include "base_particles.hpp"

namespace SPH
{
    //=================================================================================================//
SegregationLocalIsotropicDiffusion::SegregationLocalIsotropicDiffusion(const std::string &diffusion_species_name,
                                                 const std::string &gradient_species_name,
                                                 Real diff_background, Real diff_max)
    : IsotropicDiffusion(diffusion_species_name, gradient_species_name, diff_background),
      diff_max_(diff_max), local_diffusivity_(nullptr) {}
//=================================================================================================//
SegregationLocalIsotropicDiffusion::SegregationLocalIsotropicDiffusion(const std::string &species_name,
                                                 Real diff_background, Real diff_max)
    : SegregationLocalIsotropicDiffusion(species_name, species_name, diff_background, diff_max) {}
//=================================================================================================//
SegregationLocalIsotropicDiffusion::SegregationLocalIsotropicDiffusion(ConstructArgs<std::string, Real, Real> args)
    : SegregationLocalIsotropicDiffusion(std::get<0>(args), std::get<1>(args), std::get<2>(args)) {}
//=================================================================================================//
void SegregationLocalIsotropicDiffusion::initializeLocalParameters(BaseParticles *base_particles)
{
    local_diffusivity_ = base_particles->registerStateVariable<Real>(
        "SegregationDiffusivity", [&](size_t i) -> Real
        { return diff_cf_; });
    base_particles->addVariableToWrite<Real>("SegregationDiffusivity");
    
    local_segregation_rate_ = base_particles->registerStateVariable<Real>(
        "SegregationRate", [&](size_t i) -> Real
        { return diff_max_; });
    base_particles->addVariableToWrite<Real>("SegregationRate");
    
    segregation_test_ = base_particles->registerStateVariable<Real>(
        "SegregationTest");
    base_particles->addVariableToWrite<Real>("SegregationTest");

    /*Monitors*/
    segregation_component_C1_ = base_particles->registerStateVariable<Real>("ComponentC1");
    segregation_component_C2_ = base_particles->registerStateVariable<Real>("ComponentC2");
    diffusion_component_C3_ = base_particles->registerStateVariable<Real>("ComponentC3");
    segregation_component_C1_plus_C2 = base_particles->registerStateVariable<Real>("ComponentC1PlusC2");
    base_particles->addVariableToWrite<Real>("ComponentC1");
    base_particles->addVariableToWrite<Real>("ComponentC2");
    base_particles->addVariableToWrite<Real>("ComponentC3");
    base_particles->addVariableToWrite<Real>("ComponentC1PlusC2");
}
}  //namespace SPH
