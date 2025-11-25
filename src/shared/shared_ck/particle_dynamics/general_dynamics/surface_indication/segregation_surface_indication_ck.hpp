#ifndef SEGREGATION_SURFACE_INDICATION_CK_HPP
#define SEGREGATION_SURFACE_INDICATION_CK_HPP

#include "base_particles.hpp"
#include "segregation_surface_indication_ck.h"

namespace SPH
{
namespace fluid_dynamics
{

//=================================================================================================//
// SegFreeSurfaceIndicationCK<Base, RelationType<Parameters...>>
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class BaseRelationType>
SegFreeSurfaceIndicationCK<Base, RelationType<Parameters...>>::
    SegFreeSurfaceIndicationCK(BaseRelationType &base_relation)
    : Interaction<RelationType<Parameters...>>(base_relation),
      dv_indicator_(this->particles_->template registerStateVariableOnly<int>("Indicator")),
      dv_pos_div_(this->particles_->template registerStateVariableOnly<Real>("PositionDivergence")),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_threshold_by_dimensions_(0.75 * Dimensions),
      dv_smoothing_length_(this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength())
{
}

template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
SegFreeSurfaceIndicationCK<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                      SegFreeSurfaceIndicationCK<Base, RelationType<Parameters...>> &encloser,
                      Args &&...args)
    : Interaction<RelationType<Parameters...>>::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
      pos_div_(encloser.dv_pos_div_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      threshold_by_dimensions_(encloser.dv_threshold_by_dimensions_),
      smoothing_length_(encloser.dv_smoothing_length_)
{
}

//=================================================================================================//
// SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>
//=================================================================================================//
template <class FlowType, typename... Parameters>
SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>::
    SegFreeSurfaceIndicationCK(Relation<Inner<Parameters...>> &inner_relation)
    : SegFreeSurfaceIndicationCK<Base, Inner<Parameters...>>(inner_relation),
      dv_previous_surface_indicator_(
          this->particles_->template registerStateVariableOnly<int>("PreviousSurfaceIndicator"))
{
}

template <class FlowType, typename... Parameters>
template <class ExecutionPolicy>
SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                      SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>> &encloser)
    : SegFreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy))
{
}

template <class FlowType, typename... Parameters>
void SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        pos_div -= this->dW_ij(index_i, index_j) * this->Vol_[index_j] * r_ij;
    }
    this->pos_div_[index_i] = pos_div;
}

template <class FlowType, typename... Parameters>
template <class ExecutionPolicy>
SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy,
                    SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>> &encloser)
    : SegFreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy)),
      outer_(&encloser)
{
}

template <class FlowType, typename... Parameters>
void SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    bool is_near_surface = false;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        if ((this->pos_div_[index_j] < this->threshold_by_dimensions_) &&
            (r_ij < this->smoothing_length_))
        {
            is_near_surface = true;
            break;
        }
    }

    bool is_near_previous_surface = false;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        if (this->previous_surface_indicator_[index_j] == 1)
        {
            is_near_previous_surface = true;
            break;
        }
    }

    if ((this->pos_div_[index_i] < this->threshold_by_dimensions_) &&
        !is_near_surface && !is_near_previous_surface)
    {
        this->pos_div_[index_i] = 2.0 * this->threshold_by_dimensions_;
    }

    int new_indicator = 0;  // 默认非表面 🔧

    // 判断是否是表面粒子 🔧
    if ((this->pos_div_[index_i] < this->threshold_by_dimensions_))
    {
        new_indicator = 1;
    }
    else
    {
        // 判断是否靠近表面粒子 🔧
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            const UnsignedInt index_j = this->neighbor_index_[n];
            if (this->previous_surface_indicator_[index_j] == 1)
            {
                new_indicator = 2;
                break;
            }
        }
    }

    this->indicator_[index_i] = new_indicator;
    this->previous_surface_indicator_[index_i] = new_indicator == 1 ? 1 : 0; // 🔧 保持历史标记只记录表面层
}


//=================================================================================================//
// SegFreeSurfaceIndicationCK<Contact<Parameters...>>
//=================================================================================================//
template <typename... Parameters>
SegFreeSurfaceIndicationCK<Contact<Parameters...>>::
    SegFreeSurfaceIndicationCK(Relation<Contact<Parameters...>> &contact_relation)
    : SegFreeSurfaceIndicationCK<Base, Contact<Parameters...>>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy>
SegFreeSurfaceIndicationCK<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                      SegFreeSurfaceIndicationCK<Contact<Parameters...>> &encloser,
                      size_t contact_index)
    : SegFreeSurfaceIndicationCK<Base, Contact<Parameters...>>::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
void SegFreeSurfaceIndicationCK<Contact<Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        pos_div -= this->dW_ij(index_i, index_j) * this->contact_Vol_[index_j] * r_ij;
    }
    this->pos_div_[index_i] += pos_div;
}

} // namespace fluid_dynamics
} // namespace SPH

#endif // SEGREGATION_SURFACE_INDICATION_CK_HPP
