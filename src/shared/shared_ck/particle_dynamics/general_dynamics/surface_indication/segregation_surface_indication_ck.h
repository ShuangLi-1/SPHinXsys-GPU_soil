#ifndef SEGREGATION_SURFACE_INDICATION_CK_H
#define SEGREGATION_SURFACE_INDICATION_CK_H

#include "base_general_dynamics.h"
#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"

namespace SPH
{
namespace fluid_dynamics
{

template <typename... RelationTypes>
class SegFreeSurfaceIndicationCK;

//=================================================================================================//
// Base relation version
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
class SegFreeSurfaceIndicationCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
public:
    template <class BaseRelationType>
    explicit SegFreeSurfaceIndicationCK(BaseRelationType &base_relation);
    virtual ~SegFreeSurfaceIndicationCK() {}

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
    public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                          SegFreeSurfaceIndicationCK<Base, RelationType<Parameters...>> &encloser,
                          Args &&...args);

        void interact(size_t index_i, Real dt = 0.0);

    protected:
        int *indicator_;
        Real *pos_div_;
        Real *Vol_;
        Real threshold_by_dimensions_;
        Real smoothing_length_;
    };

protected:
    DiscreteVariable<int> *dv_indicator_;
    DiscreteVariable<Real> *dv_pos_div_;
    DiscreteVariable<Real> *dv_Vol_;
    Real dv_threshold_by_dimensions_;
    Real dv_smoothing_length_;
};

//=================================================================================================//
// Inner relation version with "WithUpdate"
//=================================================================================================//
template <class FlowType, typename... Parameters>
class SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>>
    : public SegFreeSurfaceIndicationCK<Base, Inner<Parameters...>>
{
public:
    explicit SegFreeSurfaceIndicationCK(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~SegFreeSurfaceIndicationCK() {}

    class InteractKernel
        : public SegFreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
    public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                          SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>> &encloser);

        void interact(size_t index_i, Real dt = 0.0);

        int *previous_surface_indicator_;
    };

    class UpdateKernel
        : public SegFreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel
    {
    public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                        SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>> &encloser);

        void update(size_t index_i, Real dt = 0.0);

    protected:
        int *previous_surface_indicator_;
        SegFreeSurfaceIndicationCK<Inner<WithUpdate, FlowType, Parameters...>> *outer_;
    };

protected:
    DiscreteVariable<int> *dv_previous_surface_indicator_;
};

using SegFreeSurfaceIndicationInnerCK = SegFreeSurfaceIndicationCK<Inner<WithUpdate, Internal>>;

//=================================================================================================//
// Contact relation version
//=================================================================================================//
template <typename... Parameters>
class SegFreeSurfaceIndicationCK<Contact<Parameters...>>
    : public SegFreeSurfaceIndicationCK<Base, Contact<Parameters...>>
{
public:
    explicit SegFreeSurfaceIndicationCK(Relation<Contact<Parameters...>> &contact_relation);
    virtual ~SegFreeSurfaceIndicationCK() {}

    class InteractKernel
        : public SegFreeSurfaceIndicationCK<Base, Contact<Parameters...>>::InteractKernel
    {
    public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                          SegFreeSurfaceIndicationCK<Contact<Parameters...>> &encloser,
                          size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

    protected:
        Real *contact_Vol_;
    };

protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};

using SegFreeSurfaceIndicationComplexCK = SegFreeSurfaceIndicationCK<Inner<WithUpdate, Internal>, Contact<>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif // SEGREGATION_SURFACE_INDICATION_CK_H
