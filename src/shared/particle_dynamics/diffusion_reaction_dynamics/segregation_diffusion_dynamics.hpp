#ifndef SEGREGATION_DIFFUSION_DYNAMICS_HPP
#define SEGREGATION_DIFFUSION_DYNAMICS_HPP

#include "segregation_diffusion_dynamics.h"

namespace SPH
{
template <class DataDelegationType, class DiffusionType>
template <class BodyRelationType>
SegregationDiffusionRelaxation<DataDelegationType, DiffusionType>::
    SegregationDiffusionRelaxation(BodyRelationType &body_relation)
    : LocalDynamics(body_relation.getSPHBody()), DataDelegationType(body_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      indicator_(this->particles_->template getVariableDataByName<int>("Indicator")),
      seg_n_(this->particles_->template getVariableDataByName<Vecd>("SegNormalDirection")),
      seg_phi_(this->particles_->template getVariableDataByName<Real>("SegPhi"))
{
    getDiffusions();
    A_ = 1.0;
    k_ = 0.0;
    for (auto &diffusion : diffusions_)
    {
        std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
        diffusion_species_.push_back(this->particles_->template registerStateVariable<Real>(diffusion_species_name));
        this->particles_->template addEvolvingVariable<Real>(diffusion_species_name);
        this->particles_->template addVariableToWrite<Real>(diffusion_species_name);
        diffusion_dt_.push_back(this->particles_->template registerStateVariable<Real>(diffusion_species_name + "ChangeRate"));

        std::string gradient_species_name = diffusion->GradientSpeciesName();
        gradient_species_.push_back(this->particles_->template registerStateVariable<Real>(gradient_species_name));
        segregation_rate_.push_back(this->particles_->template getVariableDataByName<Real>("SegregationRate"));
        segregation_test_.push_back(this->particles_->template getVariableDataByName<Real>("SegregationTest"));
        segregation_component_C1_.push_back(this->particles_->template getVariableDataByName<Real>("ComponentC1"));
        segregation_component_C2_.push_back(this->particles_->template getVariableDataByName<Real>("ComponentC2"));
        diffusion_component_C3_.push_back(this->particles_->template getVariableDataByName<Real>("ComponentC3"));
        segregation_component_C1_plus_C2_.push_back(this->particles_->template getVariableDataByName<Real>("ComponentC1PlusC2"));
        this->particles_->template addEvolvingVariable<Real>(gradient_species_name);
        this->particles_->template addVariableToWrite<Real>(gradient_species_name);
    }
}
//=================================================================================================//
template <class DataDelegationType, class DiffusionType>
void SegregationDiffusionRelaxation<DataDelegationType, DiffusionType>::getDiffusions()
{
    AbstractDiffusion &abstract_diffusion = DynamicCast<AbstractDiffusion>(this, this->sph_body_.getBaseMaterial());
    StdVec<AbstractDiffusion *> all_diffusions = abstract_diffusion.AllDiffusions();
    for (auto &diffusion : all_diffusions)
    {
        diffusions_.push_back(DynamicCast<DiffusionType>(this, diffusion));
    }
}
//=================================================================================================//
template <class DataDelegationType, class DiffusionType>
void SegregationDiffusionRelaxation<DataDelegationType, DiffusionType>::initialization(size_t index_i, Real dt)
{
    for (size_t m = 0; m < diffusions_.size(); ++m)
    {
        diffusion_dt_[m][index_i] = 0;

        segregation_component_C1_[m][index_i] = 0.0;
        segregation_component_C2_[m][index_i] = 0.0;
        segregation_component_C1_plus_C2_[m][index_i] = 0.0;
        diffusion_component_C3_[m][index_i] = 0.0;
    }
}
//=================================================================================================//
template <class DataDelegationType, class DiffusionType>
void SegregationDiffusionRelaxation<DataDelegationType, DiffusionType>::update(size_t index_i, Real dt)
{
    for (size_t m = 0; m < diffusions_.size(); ++m)
    {
        diffusion_species_[m][index_i] += dt * diffusion_dt_[m][index_i];
    }
}
//=================================================================================================//
template <class KernelGradientType, class DiffusionType>
template <typename... Args>
SegregationDiffusionRelaxation<Inner<KernelGradientType>, DiffusionType>::
    SegregationDiffusionRelaxation(Args &&...args)
    : SegregationDiffusionRelaxation<DataDelegateInner, DiffusionType>(std::forward<Args>(args)...),
      kernel_gradient_(this->particles_) {}
//=================================================================================================//
template <class KernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Inner<KernelGradientType>, DiffusionType>::interaction(size_t index_i, Real dt)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real seg_source(0.0), convec_source(0.0), q_source(0.0);
        Real C1_sum(0.0), C2_sum(0.0);
        Real test(0.0);
        auto diffusion_m = this->diffusions_[m];
        Real *gradient_species = this->gradient_species_[m];
        Real *segregation_rate = this->segregation_rate_[m];
        Real *segregation_test = this->segregation_test_[m];
        Real *segregation_component_C1 = this->segregation_component_C1_[m];
        Real *segregation_component_C2 = this->segregation_component_C2_[m];
        Real *segregation_component_C1_plus_C2 = this->segregation_component_C1_plus_C2_[m];
        Real *diffusion_component_C3 = this->diffusion_component_C3_[m];
        Real d_species = 0.0;
        Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
            Real r_ij_ = inner_neighborhood.r_ij_[n];
            Vecd &e_ij = inner_neighborhood.e_ij_[n];

            Real seg_phi_i = this->seg_phi_[index_i];
            Real seg_phi_j = this->seg_phi_[index_j];
            Real diff_coeff_ij = diffusion_m->getInterParticleDiffusionCoeff(index_i, index_j, seg_phi_i, seg_phi_j, e_ij);
            const Vecd &grad_ijV_j = this->kernel_gradient_(index_i, index_j, dW_ijV_j, e_ij);
            Real surface_area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
            Real phi_ij = gradient_species[index_i] - gradient_species[index_j];
            //d_species += diff_coeff_ij * phi_ij * surface_area_ij;
            convec_source += diff_coeff_ij * phi_ij * surface_area_ij;

            /*Segregation source*/
            Real segregation_rate_i = segregation_rate[index_i];
            Real segregation_rate_ij = segregation_rate[index_i] - segregation_rate[index_j];
            Real concentration_i =  gradient_species[index_i];
            Real F_i = this->A_*concentration_i*(1-concentration_i)*(1-this->k_*concentration_i);
            Real FF_i = 3 * this->k_*concentration_i*concentration_i - 2*(1+this->k_)*concentration_i + 1;
            FF_i *= this->A_;

            Real segregation_rate_j = segregation_rate[index_j];
            Real concentration_j =  gradient_species[index_j];
            Real F_j = this->A_*concentration_j*(1-concentration_j)*(1-this->k_*concentration_j);
            Real FF_j = 3 * this->k_*concentration_j*concentration_j - 2*(1+this->k_)*concentration_j + 1;
            FF_j *= this->A_;

            Real term1 = segregation_rate_i * FF_i * phi_ij;
            //Real term1 = segregation_rate_i * (F_i - F_j);
            Real term2 = F_i * segregation_rate_ij;
            //Real term2 = F_i* 0.1 * phi_ij;
            Real Q = segregation_rate_i * F_i - segregation_rate_j * F_j;
            Vecd gravity_dir = Vecd::Zero();
            gravity_dir[1] = -1.0;
            seg_source +=  (term1 + term2) * grad_ijV_j.dot(gravity_dir);
            //seg_source +=  Q * grad_ijV_j.dot(gravity_dir);

            /*Free-surface source term*/
            if(this->indicator_[index_i] == 1 || this->indicator_[index_i] ==2)
            {
                Real func_c_i = - segregation_rate_i * F_i * this->seg_n_[index_i].dot(gravity_dir);
                q_source += 2 * func_c_i * this->seg_n_[index_i].dot(e_ij)*dW_ijV_j;
            }  

            C1_sum += term1 * grad_ijV_j.dot(gravity_dir);
            C2_sum += term2 * grad_ijV_j.dot(gravity_dir);
        }
        this->diffusion_dt_[m][index_i] += convec_source + seg_source + q_source;
        segregation_component_C1[index_i] += C1_sum;
        segregation_component_C2[index_i] += C2_sum;
        segregation_component_C1_plus_C2[index_i] += seg_source + q_source;
        diffusion_component_C3[index_i] += convec_source;
        //segregation_test[index_i] = convec_source / (seg_source + TinyReal);
        segregation_test[index_i] = q_source;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
template <typename... Args>
SegregationDiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>::
    SegregationDiffusionRelaxation(Args &&...args)
    : SegregationDiffusionRelaxation<DataDelegateContact, DiffusionType>(
          std::forward<Args>(args)...)
{
    contact_transfer_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        contact_kernel_gradients_.push_back(ContactKernelGradientType(this->particles_, contact_particles_k));
        contact_Vol_.push_back(contact_particles_k->template registerStateVariable<Real>("VolumetricMeasure"));

        std::string diffusion_direction = "From" + this->contact_bodies_[k]->getName();
        for (auto &diffusion : this->diffusions_)
        {
            std::string variable_name = diffusion->GradientSpeciesName() + "Transfer" + diffusion_direction;
            contact_transfer_[k].push_back(
                this->particles_->template registerStateVariable<Real>(variable_name));
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>::
    resetContactTransfer(size_t index_i)
{
    for (size_t k = 0; k < this->contact_particles_.size(); ++k)
    {
        for (size_t m = 0; m < this->diffusions_.size(); ++m)
        {
            this->contact_transfer_[k][m][index_i] = 0.0;
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>::
    accumulateDiffusionRate(size_t index_i)
{
    for (size_t k = 0; k < this->contact_particles_.size(); ++k)
    {
        for (size_t m = 0; m < this->diffusions_.size(); ++m)
        {
            this->diffusion_dt_[m][index_i] += this->contact_transfer_[k][m][index_i];
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
template <typename... Args>
SegregationDiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    SegregationDiffusionRelaxation(Args &&...args)
    : SegregationDiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>(std::forward<Args>(args)...)
{
    contact_gradient_species_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        for (auto &diffusion : this->diffusions_)
        {
            std::string gradient_species_name = diffusion->GradientSpeciesName();
            contact_gradient_species_[k].push_back(
                contact_particles_k->template registerStateVariable<Real>(gradient_species_name));
            contact_particles_k->template addVariableToWrite<Real>(gradient_species_name);
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    getDiffusionChangeRateDirichlet(size_t particle_i, size_t particle_j, Vecd &e_ij,
                                    Real surface_area_ij, const StdVec<Real *> &gradient_species_k)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real diff_coeff_ij =
            this->diffusions_[m]->getInterParticleDiffusionCoeff(particle_i, particle_i, e_ij);
        Real phi_ij = 2.0 * (this->gradient_species_[m][particle_i] - gradient_species_k[m][particle_j]);
        this->diffusion_dt_[m][particle_i] += diff_coeff_ij * phi_ij * surface_area_ij;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdVec<Real *> &gradient_species_k = this->contact_gradient_species_[k];
        Real *wall_Vol_k = this->contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij_ = contact_neighborhood.r_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j, e_ij);
            Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
            getDiffusionChangeRateDirichlet(index_i, index_j, e_ij, area_ij, gradient_species_k);
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
template <typename... Args>
SegregationDiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>::
    SegregationDiffusionRelaxation(Args &&...args)
    : SegregationDiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>(
          std::forward<Args>(args)...),
      n_(this->particles_->template getVariableDataByName<Vecd>("NormalDirection"))
{
    contact_diffusive_flux_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        contact_n_.push_back(this->contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));
        contact_seg_n_.push_back(this->contact_particles_[k]->template getVariableDataByName<Vecd>("SegNormalDirection"));

        for (auto &diffusion : this->diffusions_)
        {
            std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
            contact_diffusive_flux_[k].push_back(
                contact_particles_k->template registerStateVariable<Real>(diffusion_species_name + "Flux"));
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>::
    getDiffusionChangeRateNeumann(size_t particle_i, size_t particle_j,
                                  Real surface_area_ij_Neumann,
                                  const StdVec<Real *> &diffusive_flux_k)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        this->diffusion_dt_[m][particle_i] += surface_area_ij_Neumann * diffusive_flux_k[m][particle_j];
        // this->diffusion_component_C3_[m][particle_i] += surface_area_ij_Neumann * diffusive_flux_k[m][particle_j];
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>::
    getDiffusionChangeRateNeumannInhomo(size_t particle_i, size_t particle_j,
                                  Real surface_area_ij_Neumann)
{
    Real q_source = 0.0;
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real segregation_rate_i = this->segregation_rate_[m][particle_i];
        Real concentration_i =  this->gradient_species_[m][particle_i];
        Real F_i = this->A_*concentration_i*(1-concentration_i)*(1-this->k_*concentration_i);
        Vecd gravity_dir = Vecd::Zero();
        gravity_dir[1] = -1.0;
        q_source = segregation_rate_i * F_i * n_[particle_i].dot(gravity_dir);

        this->diffusion_dt_[m][particle_i] += surface_area_ij_Neumann * q_source;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>::
    getCsmSource(size_t particle_i, size_t particle_j, Vecd contact_seg_n_j, Vecd e_ij, Real dW_ijV_j)
{
    Real csm_source(0.0);
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real segregation_rate_i = this->segregation_rate_[m][particle_i];
        Real concentration_i =  this->gradient_species_[m][particle_i];
        Real F_i = this->A_*concentration_i*(1-concentration_i)*(1-this->k_*concentration_i);
        Vecd gravity_dir = Vecd::Zero();
        gravity_dir[1] = -1.0;
        Vecd seg_n_i = this->seg_n_[particle_i];
        Real func_c_i = - segregation_rate_i * F_i * seg_n_i.dot(gravity_dir);
        //csm_source -= func_c_i * (seg_n_i + contact_seg_n_j).dot(e_ij)*dW_ijV_j;
        csm_source -= func_c_i * (contact_seg_n_j + contact_seg_n_j).dot(e_ij)*dW_ijV_j;
        this->diffusion_dt_[m][particle_i] += csm_source;
        this->segregation_component_C1_plus_C2_[m][particle_i] += csm_source;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void SegregationDiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdVec<Real *> &diffusive_flux_k = contact_diffusive_flux_[k];
        Vecd *n_k = contact_n_[k];
        Vecd *seg_n_k = contact_seg_n_[k];
        Real *Vol_k = this->contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j, e_ij);
            Vecd n_ij = n_[index_i] - n_k[index_j];
            Real area_ij_Neumann = grad_ijV_j.dot(n_ij);
            getDiffusionChangeRateNeumann(index_i, index_j, area_ij_Neumann, diffusive_flux_k);
            //getCsmSource(index_i, index_j, seg_n_k[index_j], e_ij, dW_ijV_j); 
            getCsmSource(index_i, index_j, n_k[index_j], e_ij, dW_ijV_j); 

            //getDiffusionChangeRateNeumannInhomo(index_i, index_j, area_ij_Neumann);
        }
    }
}
//=================================================================================================//
template <class SegregationDiffusionRelaxationType>
template <typename... Args>
SegregationRungeKuttaStep<SegregationDiffusionRelaxationType>::SegregationRungeKuttaStep(Args &&...args)
    : SegregationDiffusionRelaxationType(std::forward<Args>(args)...)
{
    for (auto &diffusion : this->diffusions_)
    {
        std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
        diffusion_species_s_.push_back(
            this->particles_->template registerStateVariable<Real>(diffusion_species_name + "Intermediate"));
    }
}
//=================================================================================================//
template <class SegregationDiffusionRelaxationType>
template <typename... Args>
SegregationFirstStageRK2<SegregationDiffusionRelaxationType>::SegregationFirstStageRK2(Args &&...args)
    : RungeKuttaStep<SegregationDiffusionRelaxationType>(std::forward<Args>(args)...) {}
//=================================================================================================//
template <class SegregationDiffusionRelaxationType>
void SegregationFirstStageRK2<SegregationDiffusionRelaxationType>::initialization(size_t index_i, Real dt)
{
    SegregationDiffusionRelaxationType::initialization(index_i, dt);

    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        this->diffusion_species_s_[m][index_i] = this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class SegregationDiffusionRelaxationType>
template <typename... Args>
SegregationSecondStageRK2<SegregationDiffusionRelaxationType>::SegregationSecondStageRK2(Args &&...args)
    : RungeKuttaStep<SegregationDiffusionRelaxationType>(std::forward<Args>(args)...) {}
//=================================================================================================//
template <class SegregationDiffusionRelaxationType>
void SegregationSecondStageRK2<SegregationDiffusionRelaxationType>::update(size_t index_i, Real dt)
{
    SegregationDiffusionRelaxationType::update(index_i, dt);
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real concentration_i = 0.5 * this->diffusion_species_s_[m][index_i] + 0.5 * this->diffusion_species_[m][index_i];
        concentration_i = SMIN(concentration_i, Real(1.0));
        concentration_i = SMAX(concentration_i, Real(0.0));
        this->diffusion_species_[m][index_i] = concentration_i;
    }
}
//=================================================================================================//
template <class SegregationDiffusionRelaxationType>
template <typename FirstArg, typename... OtherArgs>
SegregationDiffusionRelaxationRK2<SegregationDiffusionRelaxationType>::
    SegregationDiffusionRelaxationRK2(FirstArg &first_arg, OtherArgs &&...other_args)
    : BaseDynamics<void>(),
      rk2_1st_stage_(first_arg, std::forward<OtherArgs>(other_args)...),
      rk2_2nd_stage_(first_arg, std::forward<OtherArgs>(other_args)...) {}
//=================================================================================================//
template <class SegregationDiffusionRelaxationType>
void SegregationDiffusionRelaxationRK2<SegregationDiffusionRelaxationType>::exec(Real dt)
{
    rk2_1st_stage_.exec(dt);
    rk2_2nd_stage_.exec(dt);
}
//=================================================================================================//
} // namespace SPH
#endif // SEGREGATION_DIFFUSION_DYNAMICS_HPP