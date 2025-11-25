#ifndef SEGREGATION_DIFFUSION_CK_H
#define SEGREGATION_DIFFUSION_CK_H
#include "diffusion_dynamics_ck.h"
#include "diffusion_dynamics_ck.hpp"

namespace SPH
{
template <typename... InteractionTypes>
class SegregationDiffusionRelaxationCK;

template <class DiffusionType, class BaseInteractionType>
class SegregationDiffusionRelaxationCK<DiffusionType, BaseInteractionType>
    : public BaseInteractionType
{
    StdVec<DiffusionType *> obtainConcreteDiffusions(AbstractDiffusion &abstract_diffusion);
    StdVec<std::string> obtainDiffusionSpeciesNames(StdVec<DiffusionType *> &diffusions);
    StdVec<std::string> obtainGradientSpeciesNames(StdVec<DiffusionType *> &diffusions);

  public:
    template <class DynamicsIdentifier>
    SegregationDiffusionRelaxationCK(DynamicsIdentifier &identifier, AbstractDiffusion *abstract_diffusion);
    template <class DynamicsIdentifier>
    explicit SegregationDiffusionRelaxationCK(DynamicsIdentifier &identifier);
    template <typename BodyRelationType, typename FirstArg>
    SegregationDiffusionRelaxationCK(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : SegregationDiffusionRelaxationCK(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~SegregationDiffusionRelaxationCK() {};
    StdVec<DiffusionType *> &getDiffusions() { return diffusions_; };
    StdVec<std::string> &getDiffusionSpeciesNames() { return diffusion_species_names_; };
    StdVec<std::string> &getGradientSpeciesNames() { return gradient_species_names_; };
    DiscreteVariableArray<Real> &dvGradientSpeciesArray() { return dv_gradient_species_array_; };

    class InteractKernel : public BaseInteractionType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        DataArray<Real> *diffusion_species_, *gradient_species_, *diffusion_dt_;
        DataArray<Real> *segregation_rate_;
        DataArray<Vecd> *seg_n_;
        DataArray<Real> *segregation_phi_, *segregation_test_;
        UnsignedInt number_of_species_;
    };

  protected:
    StdVec<DiffusionType *> diffusions_;
    StdVec<std::string> diffusion_species_names_;
    StdVec<std::string> gradient_species_names_;
    DiscreteVariableArray<Real> dv_diffusion_species_array_;
    DiscreteVariableArray<Real> dv_gradient_species_array_;
    DiscreteVariableArray<Real> dv_diffusion_dt_array_;

    /*Segregation parameters*/
    StdVec<std::string> segregation_rate_names_ = {"SegregationRate"};
    DiscreteVariableArray<Real> dv_segregation_rate_array_;
    StdVec<std::string> segregation_n_names_ = {"SegNormalDirection"};
    DiscreteVariableArray<Vecd> dv_segregation_n_array_;
    StdVec<std::string> segregation_phi_names_ = {"SegregationPhi"};
    DiscreteVariableArray<Real> dv_segregation_phi_array_;
    StdVec<std::string> segregation_test_names_ = {"SegregationTest"};
    DiscreteVariableArray<Real> dv_segregation_test_array_;
};

template <class DiffusionType, class KernelGradientType, class... Parameters>
class SegregationDiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelGradientType, Parameters...>>
    : public SegregationDiffusionRelaxationCK<DiffusionType, Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = SegregationDiffusionRelaxationCK<DiffusionType, Interaction<Inner<Parameters...>>>;
    using GradientKernel = typename KernelGradientType::ComputingKernel;
    using InterParticleDiffusionCoeff = typename DiffusionType::InterParticleDiffusionCoeff;

  public:
    template <typename... Args>
    SegregationDiffusionRelaxationCK(Args &&...args);
    virtual ~SegregationDiffusionRelaxationCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        GradientKernel gradient_;
        InterParticleDiffusionCoeff *inter_particle_diffusion_coeff_;
        Real *Vol_;
        int *indicator_;
        Real smoothing_length_sq_;
    };

  protected:
    KernelGradientType kernel_gradient_;
    ConstantArray<DiffusionType, InterParticleDiffusionCoeff> ca_inter_particle_diffusion_coeff_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<int> *dv_indicator_;
    Real smoothing_length_sq_;
};

template <class DiffusionType, template <typename...> class BoundaryType, class KernelGradientType>
class SegregationDiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelGradientType>>
    : public SegregationDiffusionRelaxationCK<DiffusionType, Interaction<Contact<>>>
{
    UniquePtrsKeeper<DiscreteVariableArray<Real>> contact_transfer_array_ptrs_keeper_;
    UniquePtrsKeeper<KernelGradientType> kernel_gradient_ptrs_keeper_;
    UniquePtrsKeeper<BoundaryType<DiffusionType>> boundary_ptrs_keeper_;
    using BaseInteraction = SegregationDiffusionRelaxationCK<DiffusionType, Interaction<Contact<>>>;
    using GradientKernel = typename KernelGradientType::ComputingKernel;
    using BoundaryKernel = typename BoundaryType<DiffusionType>::ComputingKernel;

  public:
    template <typename... Args>
    explicit SegregationDiffusionRelaxationCK(Args &&...args);
    virtual ~SegregationDiffusionRelaxationCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                       UnsignedInt contact_index);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataArray<Real> *contact_transfer_;
        GradientKernel gradient_;
        BoundaryKernel boundary_flux_;

        Vecd *contact_seg_n_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
    StdVec<DiscreteVariableArray<Real> *> contact_dv_transfer_array_;
    StdVec<KernelGradientType *> contact_kernel_gradient_method_;
    StdVec<BoundaryType<DiffusionType> *> contact_boundary_method_;

    StdVec<DiscreteVariable<Vecd> *> dv_contact_seg_n_;
};

template <template <typename...> class RelationType, class... InteractionParameters>
class SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
    : public SegregationDiffusionRelaxationCK<RelationType<InteractionOnly, InteractionParameters...>>
{
    using BaseDynamicsType = SegregationDiffusionRelaxationCK<RelationType<InteractionOnly, InteractionParameters...>>;

  public:
    template <typename... Args>
    SegregationDiffusionRelaxationCK(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~SegregationDiffusionRelaxationCK() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(UnsignedInt index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_dt_;
        UnsignedInt number_of_species_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(UnsignedInt index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_, *diffusion_dt_;
        UnsignedInt number_of_species_;
    };
};

template <template <typename...> class RelationType, class... InteractionParameters>
class SegregationDiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>
    : public SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
{
    using BaseDynamicsType = SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>;

  public:
    template <typename... Args>
    SegregationDiffusionRelaxationCK(Args &&...args);
    virtual ~SegregationDiffusionRelaxationCK() {};

    class InitializeKernel : public BaseDynamicsType::InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(UnsignedInt index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_;
        DataArray<Real> *diffusion_species_s_;
    };

  protected:
    DiscreteVariableArray<Real> dv_diffusion_species_array_s_;
};

template <template <typename...> class RelationType, class... InteractionParameters>
class SegregationDiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>
    : public SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
{
    using BaseDynamicsType = SegregationDiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>;

  public:
    template <typename... Args>
    SegregationDiffusionRelaxationCK(Args &&...args);
    
    virtual ~SegregationDiffusionRelaxationCK() {};

    class UpdateKernel : public BaseDynamicsType::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(UnsignedInt index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_s_;
    };

  protected:
    DiscreteVariableArray<Real> dv_diffusion_species_array_s_;
};
} //namespace SPH

#endif //SEGREGATION_DIFFUSION_CK_H