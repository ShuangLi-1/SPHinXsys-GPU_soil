#ifndef SEGREGATION_DIFFUSION_DYNAMICS_H
#define SEGREGATION_DIFFUSION_DYNAMICS_H

#include "general_diffusion_reaction_dynamics.h"
#include "diffusion_dynamics.h"

namespace SPH
{
template <typename... InteractionTypes>
class SegregationDiffusionRelaxation;

template <class DataDelegationType, class DiffusionType>
class SegregationDiffusionRelaxation<DataDelegationType, DiffusionType>
    : public LocalDynamics,
      public DataDelegationType
{
  protected:
    Real *Vol_;
    int *indicator_;
    Vecd *seg_n_;
    Real *seg_phi_;
    StdVec<DiffusionType *> diffusions_;
    StdVec<Real *> diffusion_species_;
    StdVec<Real *> gradient_species_;
    StdVec<Real *> segregation_rate_;
    StdVec<Real *> segregation_test_;
    StdVec<Real *> segregation_component_C1_;
    StdVec<Real *> segregation_component_C2_;
    StdVec<Real *> diffusion_component_C3_;
    StdVec<Real *> segregation_component_C1_plus_C2_;
    StdVec<Real *> diffusion_dt_;
    /*Scalars*/
    Real A_, k_;



  public:
    template <class BodyRelationType>
    explicit SegregationDiffusionRelaxation(BodyRelationType &body_relation);
    void initialization(size_t index_i, Real dt = 0.0); // for contact diffusion integrated independently.
    void update(size_t index_i, Real dt = 0.0);

  private:
    void getDiffusions();
};

/**
 * @class SegregationDiffusionRelaxationInner
 * @brief Compute the diffusion relaxation process of all species
 */
template <class KernelGradientType, class DiffusionType>
class SegregationDiffusionRelaxation<Inner<KernelGradientType>, DiffusionType>
    : public SegregationDiffusionRelaxation<DataDelegateInner, DiffusionType>
{
  protected:
    KernelGradientType kernel_gradient_;

  public:
    template <typename... Args>
    explicit SegregationDiffusionRelaxation(Args &&...args);

    virtual ~SegregationDiffusionRelaxation() {};
    inline void interaction(size_t index_i, Real dt = 0.0);
};

template <class ContactKernelGradientType, class DiffusionType>
class SegregationDiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
    : public SegregationDiffusionRelaxation<DataDelegateContact, DiffusionType>
{
  protected:
    StdVec<ContactKernelGradientType> contact_kernel_gradients_;
    StdVec<Real *> contact_Vol_;
    StdVec<StdVec<Real *>> contact_transfer_;

    void resetContactTransfer(size_t index_i);
    void accumulateDiffusionRate(size_t index_i);

  public:
    template <typename... Args>
    explicit SegregationDiffusionRelaxation(Args &&...args);
    virtual ~SegregationDiffusionRelaxation() {};
};

template <class ContactKernelGradientType, class DiffusionType>
class SegregationDiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>
    : public SegregationDiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
{

  protected:
    StdVec<StdVec<Real *>> contact_gradient_species_;
    void getDiffusionChangeRateDirichlet(
        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
        const StdVec<Real *> &gradient_species_k);

  public:
    template <typename... Args>
    explicit SegregationDiffusionRelaxation(Args &&...args);
    virtual ~SegregationDiffusionRelaxation() {};
    inline void interaction(size_t index_i, Real dt = 0.0);
};

template <class ContactKernelGradientType, class DiffusionType>
class SegregationDiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>
    : public SegregationDiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
{
    Vecd *n_;
    StdVec<StdVec<Real *>> contact_diffusive_flux_;
    StdVec<Vecd *> contact_n_;
    StdVec<Vecd *> contact_seg_n_;

  protected:
    void getDiffusionChangeRateNeumann(size_t particle_i, size_t particle_j,
                                       Real surface_area_ij_Neumann,
                                       const StdVec<Real *> &diffusive_flux_k);
    void getDiffusionChangeRateNeumannInhomo(size_t particle_i, size_t particle_j,
                                       Real surface_area_ij_Neumann);
    void getCsmSource(size_t particle_i, size_t particle_j, Vecd contact_seg_n_j, Vecd e_ij, Real dW_ijV_j);

  public:
    template <typename... Args>
    explicit SegregationDiffusionRelaxation(Args &&...args);
    virtual ~SegregationDiffusionRelaxation() {};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class SegregationRungeKuttaStep
 * @brief A general step for runge-kutta integration scheme.
 * @details Am intermediate state for species is introduced here
 * to achieve multi-step integration.
 */
template <class SegregationDiffusionRelaxationType>
class SegregationRungeKuttaStep : public SegregationDiffusionRelaxationType
{
  protected:
    StdVec<Real *> diffusion_species_s_;

  public:
    template <typename... Args>
    SegregationRungeKuttaStep(Args &&...args);

    virtual ~SegregationRungeKuttaStep() {};
};
/**
 * @class SegregationFirstStageRK2
 * @brief The first stage of a 2nd-order runge-kutta integration scheme.
 * A intermediate state for species is introduced here to achieve multi-step integration.
 */
template <class SegregationDiffusionRelaxationType>
class SegregationFirstStageRK2 : public RungeKuttaStep<SegregationDiffusionRelaxationType>
{
  public:
    template <typename... Args>
    SegregationFirstStageRK2(Args &&...args);

    virtual ~SegregationFirstStageRK2() {};
    void initialization(size_t index_i, Real dt = 0.0);
};

/**
 * @class SegregationSecondStageRK2
 * @brief The second stage of the 2nd-order Runge-Kutta scheme.
 */
template <class SegregationDiffusionRelaxationType>
class SegregationSecondStageRK2 : public RungeKuttaStep<SegregationDiffusionRelaxationType>
{
  public:
    template <typename... Args>
    SegregationSecondStageRK2(Args &&...args);

    virtual ~SegregationSecondStageRK2() {};
    void update(size_t index_i, Real dt = 0.0);
};
/**
 * @class SegregationDiffusionRelaxationRK2
 * @brief The 2nd-order runge-kutta integration scheme.
 */
template <class SegregationDiffusionRelaxationType>
class SegregationDiffusionRelaxationRK2 : public BaseDynamics<void>
{
  protected:
    Dynamics1Level<SegregationFirstStageRK2<SegregationDiffusionRelaxationType>> rk2_1st_stage_;
    Dynamics1Level<SegregationSecondStageRK2<SegregationDiffusionRelaxationType>> rk2_2nd_stage_;

  public:
    template <typename FirstArg, typename... OtherArgs>
    explicit SegregationDiffusionRelaxationRK2(FirstArg &first_arg, OtherArgs &&...other_args);

    virtual ~SegregationDiffusionRelaxationRK2() {};

    virtual void exec(Real dt = 0.0) override;
};

/*Complete*/
template <class DiffusionType, class KernelGradientType, class ContactKernelGradientType,
          template <typename... Parameters> typename... ContactInteractionTypes>
class SegregationDiffusionBodyRelaxationComplex
    : public SegregationDiffusionRelaxationRK2<
          ComplexInteraction<SegregationDiffusionRelaxation<
                                 Inner<KernelGradientType>, ContactInteractionTypes<ContactKernelGradientType>...>,
                             DiffusionType>>
{
  public:
    template <typename FirstArg, typename... OtherArgs>
    explicit SegregationDiffusionBodyRelaxationComplex(FirstArg &&first_arg, OtherArgs &&...other_args)
        : SegregationDiffusionRelaxationRK2<
              ComplexInteraction<SegregationDiffusionRelaxation<
                                     Inner<KernelGradientType>, ContactInteractionTypes<ContactKernelGradientType>...>,
                                 DiffusionType>>(first_arg, std::forward<OtherArgs>(other_args)...){};
    virtual ~SegregationDiffusionBodyRelaxationComplex() {};
};
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_H
