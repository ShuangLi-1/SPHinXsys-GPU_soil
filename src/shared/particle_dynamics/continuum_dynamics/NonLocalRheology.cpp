#include "NonLocalRheology.hpp"
namespace SPH
{
namespace continuum_dynamics
{
//=============================================================================================//
UpdateNonLocalParameters::UpdateNonLocalParameters(SPHBody &sph_body, Real scale_factor)
    : LocalDynamics(sph_body),
      plastic_continuum_(DynamicCast<PlasticContinuum>(this, sph_body_.getBaseMaterial())),
      concentration_(particles_->getVariableDataByName<Real>("Concentration")),
      test_(particles_->getVariableDataByName<Real>("Test")),
      gama_(particles_->getVariableDataByName<Real>("ShearRate")),
      segregation_rate_(particles_->getVariableDataByName<Real>("SegregationRate")),
      diffusion_rate_(particles_->getVariableDataByName<Real>("SegregationDiffusivity")),
      inertial_num_(particles_->getVariableDataByName<Real>("InertialNumber")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      alpha_each_(particles_->getVariableDataByName<Real>("AlphaPhi")),
      Kc_each_(particles_->getVariableDataByName<Real>("KC")),
      friction_(particles_->getVariableDataByName<Real>("Friction")),
      fluidity_(particles_->getVariableDataByName<Real>("Fluidity")),
      local_fluidity_rate_(particles_->getVariableDataByName<Real>("LocalFluidityRate")),
      nonlocal_fluidity_rate_(particles_->getVariableDataByName<Real>("NonLocalFluidityRate")),
      scale_factor_(scale_factor),x_d_(0.0607), x_v_(0.21)
      {
        d_min_ = plastic_continuum_.getMinDiameter();
        d_max_ = plastic_continuum_.getMaxDiameter();

        A_ = 0.108;
        beta_ = 0.3744;
        R_ = d_max_ / d_min_;
        Fai_ = 0.6;
        C_ = 0.2712;
        epsilon_ = 2.0957;
        gravity_ = 9.8;

        t0_ = 1.0e-4;
        rho0_ = plastic_continuum_.ReferenceDensity();
        cohesion_ = plastic_continuum_.getCohesion();
        miu_s_ = plastic_continuum_.getStaticViscosity();
        miu_d_ = plastic_continuum_.getDynamicViscosity();
        I0_ = plastic_continuum_.getMaterialConstant();
      }
//=============================================================================================//
void UpdateNonLocalParameters::update(size_t index_i, Real dt)
{
    Real concentration_i = concentration_[index_i];
    Real d_avg = concentration_i * d_min_ + (1-concentration_i) * d_max_;
    test_[index_i] = d_avg;
    Real gama_i = gama_[index_i] * scale_factor_;
    Real rho_i = rho_[index_i];
    Real p_i = p_[index_i];
    Real miu_i = friction_[index_i];
    Real fluidity_i = fluidity_[index_i];
    Real fluidity_new(0.0), miu_new(0.0);
    Real local_fluidity_rate_i(0.0);
    Real nonlocal_fluidity_rate_i = nonlocal_fluidity_rate_[index_i];
    Real total_fluidity_rate(0.0);
    if (p_i < -1.0)
    {
        miu_new = TinyReal;
        fluidity_new = TinyReal;
        // miu_new = miu_d_;
        // fluidity_new = SMAX(gama_i / miu_new, TinyReal);
    }
    if (p_i >= -1.0 && p_i <= 1.0)
    {
        miu_new = miu_d_;
        fluidity_new = SMAX(gama_i / miu_new, TinyReal);
    }
    if (p_i > 1.0)
    {
        local_fluidity_rate_i = -(miu_s_ - miu_i) * fluidity_i 
        - (miu_d_ - miu_s_) * sqrt(rho0_ * d_avg * d_avg /(p_i + TinyReal)) *miu_i*fluidity_i*fluidity_i/ I0_;
        total_fluidity_rate = local_fluidity_rate_i + nonlocal_fluidity_rate_i;
        //total_fluidity_rate = local_fluidity_rate_i;

        local_fluidity_rate_[index_i] = local_fluidity_rate_i + nonlocal_fluidity_rate_i;
        fluidity_new = SMAX(TinyReal, fluidity_i + total_fluidity_rate * dt / t0_);
        miu_new = gama_i/fluidity_new;
        miu_new = clamp(miu_new, miu_s_, miu_d_);
    }
    
    fluidity_[index_i] = fluidity_new;
    friction_[index_i] = miu_new;

    //alpha_each_[index_i] = plastic_continuum_.getDPConstantsA(35.0 * Pi / 180);
    alpha_each_[index_i] = plastic_continuum_.getDPConstantsA_WithMiu(miu_new);
    //Kc_each_[index_i] = plastic_continuum_.getDPConstantsK(500,35.0 * Pi / 180);
    Kc_each_[index_i] = plastic_continuum_.getDPConstantsK_WithMiu(miu_new);

}
} // namespace continuum_dynamics
} // namespace SPH