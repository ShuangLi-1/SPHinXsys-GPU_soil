#ifndef DIFFUSION_SEGREGATION_REACTION_H
#define DIFFUSION_SEGREGATION_REACTION_H

#include "diffusion_reaction.h"


namespace SPH
{
/**
 * @class SegregationLocalIsotropicDiffusion
 * @brief diffusion coefficient D in Zhu et al(2023).
 * TODO: The difference between algebraic and geometric average should be identified.
 */
class SegregationLocalIsotropicDiffusion : public IsotropicDiffusion
{
  protected:
    Real diff_max_; /**< maximum diffusion coefficient. */
    Real *local_diffusivity_, *local_segregation_rate_, *segregation_test_;
    Real *segregation_component_C1_, *segregation_component_C2_, *diffusion_component_C3_, *segregation_component_C1_plus_C2;

  public:
    SegregationLocalIsotropicDiffusion(const std::string &diffusion_species_name,
                            const std::string &gradient_species_name,
                            Real diff_background, Real diff_max);
    SegregationLocalIsotropicDiffusion(const std::string &species_name, Real diff_background, Real diff_max);
    explicit SegregationLocalIsotropicDiffusion(ConstructArgs<std::string, Real, Real> args);
    virtual ~SegregationLocalIsotropicDiffusion() {};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Real getReferenceDiffusivity() override { return diff_max_; };
    virtual Real getDiffusionCoeffWithBoundary(size_t index_i) override { return local_diffusivity_[index_i]; };
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        return 0.5 * (local_diffusivity_[index_i] + local_diffusivity_[index_j]);
    };
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, Real seg_phi_i, Real seg_phi_j,const Vecd &e_ij)
    {
            return 0.5*(local_diffusivity_[index_i]*seg_phi_i + local_diffusivity_[index_j]*seg_phi_j);
    };

    virtual Real getSegregationRate(size_t index_i)
    {
      return local_segregation_rate_[index_i];
    }

    virtual Real getInterParticleSegregationRate(size_t index_i, size_t index_j)
    {
      return local_segregation_rate_[index_i] - local_segregation_rate_[index_j];
    }


    class InterParticleDiffusionCoeff
    {
        Real diff_cf_;

      public:
        InterParticleDiffusionCoeff() : diff_cf_(0) {};
        InterParticleDiffusionCoeff(SegregationLocalIsotropicDiffusion &encloser)
            : diff_cf_(encloser.diff_cf_) {};
        template <class ExecutionPolicy>
        InterParticleDiffusionCoeff(const ExecutionPolicy &ex_policy, SegregationLocalIsotropicDiffusion &encloser)
            : InterParticleDiffusionCoeff(encloser){};
        Real operator()(size_t index_i, size_t index_j, Real seg_phi_i, Real seg_phi_j,const Vecd &e_ij)
        {
            return 0.5*(diff_cf_*seg_phi_i + diff_cf_*seg_phi_j);
        };
        Real operator()(size_t index_i, size_t index_j,const Vecd &e_ij)
        {
            return 0.5*(diff_cf_ + diff_cf_);
        };
    };
};
}  //namespace SPH

#endif //DIFFUSION_SEGREGATION_REACTION_H