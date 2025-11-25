#ifndef CONTINUUM_SEGREGATION_PARAMETERS_HPP
#define CONTINUUM_SEGREGATION_PARAMETERS_HPP

#include "continuum_segregation_parameters.h"
#include "base_particles.hpp"
namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
template <class BaseInteractionType>
template <class DynamicsIdentifier>
SegregationParametersUpdateBase<BaseInteractionType>::SegregationParametersUpdateBase(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
    material_(DynamicCast<BaseMaterial>(this, this->sph_body_.getBaseMaterial())),
    dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
    dv_seg_phi_(this->particles_->template registerStateVariableOnly<Real>("SegregationPhi", Real(1.0))),
    dv_seg_n_(this->particles_->template registerStateVariableOnly<Vecd>("SegNormalDirection")),
    dv_indicator_(this->particles_->template getVariableByName<int>("Indicator"))
    {
        this->particles_->template addEvolvingVariable<Vecd>("SegNormalDirection");
        this->particles_->template addEvolvingVariable<Real>("SegregationPhi");
    }
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
UpdateSegNormalDirection<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateSegNormalDirection(Relation<Inner<Parameters...>> &inner_relation)
    : SegregationParametersUpdateBase<Interaction<Inner<Parameters...>>>(inner_relation),
      correction_(this->particles_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
UpdateSegNormalDirection<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : seg_n_(encloser.dv_seg_n_->DelegatedData(ex_policy)),
      seg_phi_(encloser.dv_seg_phi_->DelegatedData(ex_policy))
    {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void UpdateSegNormalDirection<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    seg_n_[index_i] = Vecd::Zero();
    seg_phi_[index_i] = 0.0;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
UpdateSegNormalDirection<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(encloser.correction_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      seg_n_(encloser.dv_seg_n_->DelegatedData(ex_policy)),
      seg_phi_(encloser.dv_seg_phi_->DelegatedData(ex_policy)),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
      W0_(this->kernel_.W(ZeroData<Vecd>::value))
      {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void UpdateSegNormalDirection<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real phi_sum_wvj(W0_ * Vol_[index_i]);
    //Real phi_sum_wvj(0.0);
    Vecd sum_dwvj_ = Vecd::Zero();
    int indicator_i = indicator_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd nablaW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j] * this->e_ij(index_i, index_j);
        Real W_ijV_j = this->W_ij(index_i, index_j) * Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        int indicator_j = indicator_[index_j];
        /*Degeneration coeff*/
        phi_sum_wvj += W_ijV_j;

        /*Segregation normal direction*/
        if(indicator_i == 1)
            sum_dwvj_ += dW_ijV_j * e_ij;
        if(indicator_i == 2 && indicator_j ==1)
            sum_dwvj_ -= dW_ijV_j * e_ij;  
    }
    seg_phi_[index_i] += phi_sum_wvj;
    seg_n_[index_i] += sum_dwvj_;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
UpdateSegNormalDirection<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : seg_n_(encloser.dv_seg_n_->DelegatedData(ex_policy)),
      seg_phi_(encloser.dv_seg_phi_->DelegatedData(ex_policy)),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)){}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void UpdateSegNormalDirection<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    seg_n_[index_i] = seg_n_[index_i] / (seg_n_[index_i].norm()+TinyReal);
    if(indicator_[index_i] == 1 || indicator_[index_i] == 2)
        seg_phi_[index_i] = 2.0 * seg_phi_[index_i] - 1.0;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
UpdateSegNormalDirection<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateSegNormalDirection(Relation<Contact<Parameters...>> &wall_contact_relation)
    : SegregationParametersUpdateBase<Interaction<Contact<Wall, Parameters...>>>(wall_contact_relation),
      correction_(this->particles_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
UpdateSegNormalDirection<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(encloser.correction_),
      wall_Vol_(encloser.dv_wall_Vol_[contact_index]->DelegatedData(ex_policy)),
      seg_n_(encloser.dv_seg_n_->DelegatedData(ex_policy)),
      seg_phi_(encloser.dv_seg_phi_->DelegatedData(ex_policy)),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void UpdateSegNormalDirection<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd contact_seg_n = Vecd::Zero();
    Real phi_WijVj = seg_phi_[index_i];
    Real sum_WijVj = seg_phi_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * wall_Vol_[index_j];
        Real W_ijV_j = this->W_ij(index_i, index_j) * wall_Vol_[index_j];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        /*Normal Direction*/
        contact_seg_n -= dW_ijV_j * e_ij;
        /*Second-order correction*/
        phi_WijVj -= W_ijV_j;
        sum_WijVj += W_ijV_j;
        // if(this->indicator_[index_i] == 0)
        //     this->indicator_[index_i] = 3;
    }
    seg_n_[index_i] += contact_seg_n;
    /*Only for particles near wall*/
    if(this->indicator_[index_i] == 3)
    {
        Real seg_phi_i = phi_WijVj/(TinyReal + sum_WijVj);
        seg_phi_[index_i] = seg_phi_i;
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
UpdateSegNormalDirection<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateSegNormalDirection(Relation<Contact<Parameters...>> &soil_contact_relation)
    : SegregationParametersUpdateBase<Interaction<Contact<Soil, Parameters...>>>(soil_contact_relation),
      correction_(this->particles_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
UpdateSegNormalDirection<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(encloser.correction_),
      soil_Vol_(encloser.dv_soil_Vol_[contact_index]->DelegatedData(ex_policy)),
      seg_n_(encloser.dv_seg_n_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void UpdateSegNormalDirection<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd contact_seg_n = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * soil_Vol_[index_j];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        contact_seg_n += dW_ijV_j * e_ij;
    }
    seg_n_[index_i] += contact_seg_n;
}

} // namespace continuum_dynamics
} // namespace SPH
#endif //CONTINUUM_SEGREGATION_PARAMETERS_HPP