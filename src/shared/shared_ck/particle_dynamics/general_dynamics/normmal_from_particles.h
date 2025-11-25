
#ifndef NORMAL_FROM_PARTICLES_H
#define NORMAL_FROM_PARTICLES_H

#include "interaction_ck.hpp"

namespace SPH
{
template <class BaseInteractionType>
class BaseNormalFromParticles : public BaseInteractionType
{

  public:
    template <class DynamicsIdentifier>
    explicit BaseNormalFromParticles(DynamicsIdentifier &identifier);
    virtual ~BaseNormalFromParticles() {};

  protected:
    Shape *initial_shape_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_n_, *dv_n0_;
    DiscreteVariable<Real> *dv_phi_, *dv_phi0_, *dv_Vol_;
};
}  //namespace SPH
#endif //NORMAL_FROM_PARTICLES_H
