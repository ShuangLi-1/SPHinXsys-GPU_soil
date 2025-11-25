#include "segregation_surface_indication.hpp"

namespace SPH
{
//=================================================================================================//
SegregationFreeSurfaceIndication<Inner<>>::
    SegregationFreeSurfaceIndication(BaseInnerRelation &inner_relation)
    : SegregationFreeSurfaceIndication<DataDelegateInner>(inner_relation),
      smoothing_length_(inner_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength()) {}
//=================================================================================================//
void SegregationFreeSurfaceIndication<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        pos_div -= inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.r_ij_[n];
    }
    pos_div_[index_i] = pos_div;
}
//=================================================================================================//
void SegregationFreeSurfaceIndication<Inner<>>::update(size_t index_i, Real dt)
{  
    int new_indicator = 0;  // 默认非表面 🔧

    // 判断是否是表面粒子 🔧
    if ((pos_div_[index_i] < threshold_by_dimensions_))
    {
        new_indicator = 1;
    }
    else if(isVeryNearFreeSurface(index_i))
    {
        new_indicator = 2;
    }

    indicator_[index_i] = new_indicator;
}
//=================================================================================================//
bool SegregationFreeSurfaceIndication<Inner<>>::isVeryNearFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        /** Two layer particles.*/
        if (pos_div_[inner_neighborhood.j_[n]] < threshold_by_dimensions_ &&
            inner_neighborhood.r_ij_[n] < 2.0 * smoothing_length_)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
//=================================================================================================//
SegregationFreeSurfaceIndication<Inner<SpatialTemporal>>::
    SegregationFreeSurfaceIndication(BaseInnerRelation &inner_relation)
    : SegregationFreeSurfaceIndication<Inner<>>(inner_relation),
      previous_surface_indicator_(particles_->registerStateVariable<int>("PreviousSurfaceIndicator", 1))
{
    particles_->addEvolvingVariable<int>("PreviousSurfaceIndicator");
}
//=================================================================================================//
void SegregationFreeSurfaceIndication<Inner<SpatialTemporal>>::interaction(size_t index_i, Real dt)
{
    SegregationFreeSurfaceIndication<Inner<>>::interaction(index_i, dt);

    if (pos_div_[index_i] < threshold_by_dimensions_ &&
        previous_surface_indicator_[index_i] != 1 &&
        !isNearPreviousFreeSurface(index_i))
        pos_div_[index_i] = 2.0 * threshold_by_dimensions_;
}
//=================================================================================================//
bool SegregationFreeSurfaceIndication<Inner<SpatialTemporal>>::isNearPreviousFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        if (previous_surface_indicator_[inner_neighborhood.j_[n]] == 1)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
//=================================================================================================//
void SegregationFreeSurfaceIndication<Inner<SpatialTemporal>>::update(size_t index_i, Real dt)
{
    SegregationFreeSurfaceIndication<Inner<>>::update(index_i, dt);

    previous_surface_indicator_[index_i] = indicator_[index_i];
}
//=================================================================================================//
void SegregationFreeSurfaceIndication<Contact<>>::interaction(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            pos_div -= contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.r_ij_[n];
        }
    }
    pos_div_[index_i] += pos_div;
}
//=================================================================================================//
} // namespace SPH
