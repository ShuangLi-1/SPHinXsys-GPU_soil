
#ifndef SEGREGATION_DIFFUSION_CK_HPP
#define SEGREGATION_DIFFUSION_CK_HPP

#include "segregation_diffusion_ck.h"

namespace SPH
{
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class DynamicsIdentifier>
SegregationDiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    SegregationDiffusionRelaxationCK(DynamicsIdentifier &identifier, AbstractDiffusion *abstract_diffusion)
    : BaseInteractionType(identifier),
      diffusions_(this->obtainConcreteDiffusions(*abstract_diffusion)),
      diffusion_species_names_(this->obtainDiffusionSpeciesNames(diffusions_)),
      gradient_species_names_(this->obtainGradientSpeciesNames(diffusions_)),
      dv_diffusion_species_array_(this->particles_->template getVariablesByName<Real>(
          diffusion_species_names_, "")),
      dv_gradient_species_array_(this->particles_->template getVariablesByName<Real>(
          gradient_species_names_, "")),
      dv_diffusion_dt_array_(this->particles_->template registerStateVariables<Real>(
          diffusion_species_names_, "ChangeRate")),
      dv_segregation_rate_array_(this->particles_->template getVariablesByName<Real>(
          segregation_rate_names_, "")),
      dv_segregation_n_array_(this->particles_->template getVariablesByName<Vecd>(
          segregation_n_names_, "")),
      dv_segregation_phi_array_(this->particles_->template getVariablesByName<Real>(
          segregation_phi_names_, "")),
      dv_segregation_test_array_(this->particles_->template getVariablesByName<Real>(
          segregation_test_names_, ""))
{
    this->particles_->template addVariableToWrite<Real>(&dv_diffusion_species_array_);
    this->particles_->template addVariableToWrite<Real>(&dv_gradient_species_array_);
    this->particles_->template addEvolvingVariable<Real>(&dv_diffusion_species_array_);
    this->particles_->template addEvolvingVariable<Real>(&dv_gradient_species_array_);

    this->particles_->template addVariableToWrite<Real>(&dv_segregation_rate_array_);
    this->particles_->template addVariableToWrite<Vecd>(&dv_segregation_n_array_);
    this->particles_->template addVariableToWrite<Real>(&dv_segregation_phi_array_);
    this->particles_->template addVariableToWrite<Real>(&dv_segregation_test_array_);
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class DynamicsIdentifier>
SegregationDiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    SegregationDiffusionRelaxationCK(DynamicsIdentifier &identifier)
    : SegregationDiffusionRelaxationCK(identifier, DynamicCast<AbstractDiffusion>(
                                            this, &identifier.getSPHBody().getBaseMaterial())) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
StdVec<DiffusionType *> SegregationDiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    obtainConcreteDiffusions(AbstractDiffusion &abstract_diffusion)
{
    StdVec<AbstractDiffusion *> all_diffusions = abstract_diffusion.AllDiffusions();
    StdVec<DiffusionType *> diffusions;
    for (auto &diffusion : all_diffusions)
    {
        diffusions.push_back(DynamicCast<DiffusionType>(this, diffusion));
    }
    return diffusions;
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
StdVec<std::string> SegregationDiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    obtainDiffusionSpeciesNames(StdVec<DiffusionType *> &diffusions)
{
    StdVec<std::string> diffusion_species_names;
    for (auto &diffusion : diffusions)
    {
        diffusion_species_names.push_back(diffusion->DiffusionSpeciesName());
    }
    return diffusion_species_names;
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
StdVec<std::string> SegregationDiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    obtainGradientSpeciesNames(StdVec<DiffusionType *> &diffusions)
{
    StdVec<std::string> gradient_species_names;
    for (auto &diffusion : diffusions)
    {
        gradient_species_names.push_back(diffusion->GradientSpeciesName());
    }
    return gradient_species_names;
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class ExecutionPolicy, class EncloserType, typename... Args>
SegregationDiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : BaseInteractionType::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      gradient_species_(encloser.dv_gradient_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      segregation_rate_(encloser.dv_segregation_rate_array_.DelegatedDataArray(ex_policy)),
      seg_n_(encloser.dv_segregation_n_array_.DelegatedDataArray(ex_policy)),
      segregation_phi_(encloser.dv_segregation_phi_array_.DelegatedDataArray(ex_policy)),
      segregation_test_(encloser.dv_segregation_test_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <class DiffusionType, class KernelGradientType, class... Parameters>
template <typename... Args>
SegregationDiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelGradientType, Parameters...>>::
    SegregationDiffusionRelaxationCK(Args &&...args)
    : BaseInteraction(std::forward<Args>(args)...),
      kernel_gradient_(this->particles_),
      ca_inter_particle_diffusion_coeff_(this->diffusions_),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_indicator_(this->particles_->template getVariableByName<int>("Indicator")),
      smoothing_length_sq_(pow(this->sph_adaptation_->ReferenceSmoothingLength(), 2)) {}
//=================================================================================================//
template <class DiffusionType, class KernelGradientType, class... Parameters>
template <class ExecutionPolicy, class EncloserType>
SegregationDiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      gradient_(ex_policy, encloser.kernel_gradient_),
      inter_particle_diffusion_coeff_(encloser.ca_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
      smoothing_length_sq_(encloser.smoothing_length_sq_) {}
//=================================================================================================//
template <class DiffusionType, class KernelGradientType, class... Parameters>
void SegregationDiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < this->number_of_species_; ++m)
    {
        Real d_species = 0.0;
        Real A_(1.0),k_(0.0);
        Real seg_source(0.0), convec_source(0.0), q_source(0.0);
        Real test(0.0);
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            Real dW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j];
            Vecd e_ij = this->e_ij(index_i, index_j);
            Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);

            /*Convection source*/
            Real surface_area_ij = 2.0 * gradient_(index_i, index_j, dW_ijV_j, e_ij).dot(vec_r_ij) /
                                   (vec_r_ij.squaredNorm() + 0.01 * this->smoothing_length_sq_);
            Real phi_ij = this->gradient_species_[m][index_i] - this->gradient_species_[m][index_j];
            Real seg_phi_i = this->segregation_phi_[m][index_i];
            Real seg_phi_j = this->segregation_phi_[m][index_j];
            convec_source += inter_particle_diffusion_coeff_[m](index_i, index_j,seg_phi_i,seg_phi_j, e_ij) * phi_ij * surface_area_ij;
            test += dW_ijV_j*phi_ij*e_ij.dot(vec_r_ij) /
                                   (vec_r_ij.squaredNorm() + 0.01 * this->smoothing_length_sq_);
            //d_species += inter_particle_diffusion_coeff_[m](index_i, index_j,1.0,1.0, e_ij) * phi_ij * surface_area_ij;
            /*Segregation source*/
            Real segregation_rate_i = this->segregation_rate_[m][index_i];
            Real segregation_rate_ij = this->segregation_rate_[m][index_i] - this->segregation_rate_[m][index_j];
            Real concentration_i =  this->gradient_species_[m][index_i];
            Real F_i = A_*concentration_i*(1-concentration_i)*(1-k_*concentration_i);
            Real FF_i = 3*k_*concentration_i*concentration_i - 2*(1+k_)*concentration_i + 1;
            FF_i *= A_;
            Real term1 = segregation_rate_i * FF_i * phi_ij;
            Real term2 = F_i * segregation_rate_ij;
            Vecd gravity_dir = Vecd::Zero();
            gravity_dir[1] = -1.0;
            seg_source +=  (term1 + term2) * gradient_(index_i, index_j, dW_ijV_j, e_ij).dot(gravity_dir);
            //d_species += (term1 + term2) * gradient_(index_i, index_j, dW_ijV_j, e_ij).dot(gravity_dir);
            
            /*Freesurface source term*/
            if(indicator_[index_i] == 1 || indicator_[index_i] ==2)
            {
                Real func_c_i = - segregation_rate_i * F_i * this->seg_n_[m][index_i].dot(gravity_dir);
                q_source += 2*func_c_i * this->seg_n_[m][index_i].dot(e_ij)*dW_ijV_j;
            }
        }
        this->diffusion_dt_[m][index_i] += convec_source ;
    }
}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelGradientType>
template <typename... Args>
SegregationDiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelGradientType>>::
    SegregationDiffusionRelaxationCK(Args &&...args)
    : BaseInteraction(std::forward<Args>(args)...)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        contact_dv_transfer_array_.push_back(
            contact_transfer_array_ptrs_keeper_.createPtr<DiscreteVariableArray<Real>>(
                this->particles_->template registerStateVariables<Real>(
                    this->diffusion_species_names_, "TransferWith" + this->sph_body_.getName())));
        contact_kernel_gradient_method_.push_back(
            kernel_gradient_ptrs_keeper_.template createPtr<KernelGradientType>(
                this->particles_, this->contact_particles_[k]));
        contact_boundary_method_.push_back(
            boundary_ptrs_keeper_.template createPtr<BoundaryType<DiffusionType>>(
                *this, this->contact_particles_[k]));

        dv_contact_seg_n_.push_back(
            this->contact_particles_[k]->template getVariableByName<Vecd>("SegNormalDirection"));
    }
}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelGradientType>
template <class ExecutionPolicy, class EncloserType>
SegregationDiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelGradientType>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_transfer_(encloser.contact_dv_transfer_array_[contact_index]->DelegatedDataArray(ex_policy)),
      gradient_(ex_policy, *encloser.contact_kernel_gradient_method_[contact_index]),
      boundary_flux_(ex_policy, *encloser.contact_boundary_method_[contact_index]),
      contact_seg_n_(encloser.dv_contact_seg_n_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelGradientType>
void SegregationDiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelGradientType>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < this->number_of_species_; ++m)
    {
        contact_transfer_[m][index_i] = 0.0;
        Real neumann_source(0.0);
        Real csm_source(0.0);
        Real A_(1.0),k_(0.0);
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            Real dW_ijV_j = this->dW_ij(index_i, index_j) * this->contact_Vol_[index_j];
            Vecd e_ij = this->e_ij(index_i, index_j);
            Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);
            
            /*Neumann Boundary*/
            Vecd surface_area_ij = 2.0 * gradient_(index_i, index_j, dW_ijV_j, e_ij);
            neumann_source += boundary_flux_(m, index_i, index_j, e_ij, vec_r_ij).dot(surface_area_ij);

            /*CSM boundary*/
            Real segregation_rate_i = this->segregation_rate_[m][index_i];
            Real concentration_i =  this->gradient_species_[m][index_i];
            Real F_i = A_*concentration_i*(1-concentration_i)*(1-k_*concentration_i);
            Vecd gravity_dir = Vecd::Zero();
            gravity_dir[1] = -1.0;
            Real func_c_i = - segregation_rate_i * F_i * this->seg_n_[m][index_i].dot(gravity_dir);
            csm_source -= func_c_i * (this->seg_n_[m][index_i] + contact_seg_n_[index_j]).dot(e_ij)*dW_ijV_j;

        }
        this->diffusion_dt_[m][index_i] += neumann_source+csm_source;

        this->segregation_test_[m][index_i] = csm_source;
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    InitializeKernel::initialize(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < number_of_species_; ++m)
    {
        diffusion_dt_[m][index_i] = 0;
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < number_of_species_; ++m)
    {
        diffusion_species_[m][index_i] += dt * diffusion_dt_[m][index_i];
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <typename... Args>
SegregationDiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    SegregationDiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_diffusion_species_array_s_(this->particles_->template registerStateVariables<Real>(
          this->diffusion_species_names_, "Intermediate")) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
SegregationDiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::InitializeKernel(ex_policy, encloser),
      diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_species_s_(encloser.dv_diffusion_species_array_s_.DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void SegregationDiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    InitializeKernel::initialize(UnsignedInt index_i, Real dt)
{
    BaseDynamicsType::InitializeKernel::initialize(index_i, dt);

    for (UnsignedInt m = 0; m < this->number_of_species_; ++m)
    {
        diffusion_species_s_[m][index_i] = diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <typename... Args>
SegregationDiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    SegregationDiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_diffusion_species_array_s_(this->particles_->template getVariablesByName<Real>(
          this->diffusion_species_names_, "Intermediate")) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
SegregationDiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::UpdateKernel(ex_policy, encloser),
      diffusion_species_s_(encloser.dv_diffusion_species_array_s_.DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void SegregationDiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    BaseDynamicsType::UpdateKernel::update(index_i, dt);
    for (UnsignedInt m = 0; m < this->number_of_species_; ++m)
    {
        Real concentration_i = 0.5 * diffusion_species_s_[m][index_i] +
                                0.5 * this->diffusion_species_[m][index_i];                               
        this->diffusion_species_[m][index_i] = concentration_i;
    }
}
} //namespace SPH

#endif //SEGREGATION_DIFFUSION_CK_HPP