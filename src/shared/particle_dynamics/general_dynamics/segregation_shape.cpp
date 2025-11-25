#include "segregation_shape.hpp"

namespace SPH
{
//=================================================================================================//
void SegregationNormalAndDegeneration<Inner<>>::initialization(size_t index_i, Real dt)
{
    seg_n_[index_i] = Vecd::Zero();
    seg_phi_[index_i] = 0.0;
}
//=================================================================================================//
void SegregationNormalAndDegeneration<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real phi_sum_wvj(W0_ * Vol_[index_i]);
    Vecd sum_dwvj_ = Vecd::Zero();
    int indicator_i = indicator_[index_i];
    Matd velocity_gradient = Matd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n]* Vol_[index_j];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        Real W_ijV_j = inner_neighborhood.W_ij_[n] * Vol_[index_j];
        Vecd e_ij =  inner_neighborhood.e_ij_[n];
        int indicator_j = indicator_[index_j];
        /*Degeneration coeff*/
        phi_sum_wvj += W_ijV_j;

        /*Segregation normal direction*/
        if(indicator_i == 1)
            sum_dwvj_ += dW_ijV_j * e_ij;
        if(indicator_i == 2 && indicator_j ==1)
            sum_dwvj_ -= dW_ijV_j * e_ij;  

        /*Velocity Gradient*/
        velocity_gradient -= (vel_[index_i] - vel_[index_j]) * dW_ijV_j * e_ij.transpose();
    }
    seg_phi_[index_i] += phi_sum_wvj;
    seg_n_[index_i] += sum_dwvj_;
    velocity_gradient_[index_i] = velocity_gradient;
}
//=================================================================================================//
void SegregationNormalAndDegeneration<Inner<>>::update(size_t index_i, Real dt)
{
    /*Segregation normal direction*/
    seg_n_[index_i] = seg_n_[index_i] / (seg_n_[index_i].norm()+TinyReal);
    /*Degeneration coefficient*/
    if(indicator_[index_i] == 1 || indicator_[index_i] == 2)
        seg_phi_[index_i] = 2.0 * seg_phi_[index_i] - 1.0;
    /*Shear rate*/
    const Real stress_dimension_ = 3.0; 
    Mat3d velocity_gradient = upgradeToMat3d(velocity_gradient_[index_i]);
    Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    Mat3d deviatoric_strain_rate = strain_rate - (1.0 / stress_dimension_) * strain_rate.trace() * Mat3d::Identity();
    Real shear_rate = sqrt(2.0 * deviatoric_strain_rate.squaredNorm());
    gama_[index_i] = shear_rate;
}
//=================================================================================================//
SegregationNormalAndDegeneration<Contact<Base>>::SegregationNormalAndDegeneration(BaseContactRelation &contact_relation)
    : SegregationNormalAndDegeneration<Base, DataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(this->contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
        contact_vel_.push_back(this->contact_particles_[k]->template getVariableDataByName<Vecd>("Velocity"));
    }
}
//=================================================================================================//
void SegregationNormalAndDegeneration<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Vecd contact_seg_n = Vecd::Zero();
    Real phi_WijVj = seg_phi_[index_i];
    Real sum_WijVj = seg_phi_[index_i];
    Matd velocity_gradient = Matd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real *contact_Vol_k = this->contact_Vol_[k];
        Vecd *cantact_vel_k = this->contact_vel_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            UnsignedInt index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * contact_Vol_k[index_j];
            Real W_ijV_j = contact_neighborhood.W_ij_[n] * contact_Vol_k[index_j];
            Real r_ij = contact_neighborhood.r_ij_[n];

            /*Normal Direction*/
            contact_seg_n -= dW_ijV_j * e_ij;
            /*Second-order correction*/
            phi_WijVj -= W_ijV_j;
            sum_WijVj += W_ijV_j;
            // if(this->indicator_[index_i] == 0)
            //     this->indicator_[index_i] = 3;
            /*Velocity Gradient*/
            Vecd vel_j_in_wall = 2.0 * cantact_vel_k[index_j] - this->vel_[index_i];
            velocity_gradient -= (this->vel_[index_i] - vel_j_in_wall) * dW_ijV_j * e_ij.transpose();
        }
    }
    seg_n_[index_i] += contact_seg_n;
    /*Only for particles near wall*/
    if(this->indicator_[index_i] == 3)
    {
        Real seg_phi_i = phi_WijVj/(TinyReal + sum_WijVj);
        seg_phi_[index_i] = seg_phi_i;
    }
    this->velocity_gradient_[index_i] += velocity_gradient;
}
//=================================================================================================//
void SegregationNormalAndDegeneration<Contact<Soil>>::interaction(size_t index_i, Real dt)
{
    Vecd contact_seg_n = Vecd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real *contact_Vol_k = this->contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            UnsignedInt index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * contact_Vol_k[index_j];

            contact_seg_n += dW_ijV_j * e_ij;
        }
    }
    seg_n_[index_i] += contact_seg_n;
}
//=================================================================================================//
//=============================================================================================//
UpdateDiffusivityAndSegregationRate::UpdateDiffusivityAndSegregationRate(SPHBody &sph_body, Real scale_factor)
    : LocalDynamics(sph_body),
      plastic_continuum_(DynamicCast<PlasticContinuum>(this, sph_body_.getBaseMaterial())),
      test_(particles_->getVariableDataByName<Real>("Test")),
      concentration_(particles_->getVariableDataByName<Real>("Concentration")),
      gama_(particles_->getVariableDataByName<Real>("ShearRate")),
      segregation_rate_(particles_->getVariableDataByName<Real>("SegregationRate")),
      diffusion_rate_(particles_->getVariableDataByName<Real>("SegregationDiffusivity")),
      inertial_num_(particles_->getVariableDataByName<Real>("InertialNumber")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      alpha_each_(particles_->getVariableDataByName<Real>("AlphaPhi")),
      Kc_each_(particles_->getVariableDataByName<Real>("KC")),
      scale_factor_(scale_factor),x_d_(0.091), x_v_(0.191)
      {
        d_min_ = plastic_continuum_.getMinDiameter() * scale_factor_;
        d_max_ = plastic_continuum_.getMaxDiameter() * scale_factor_;

        A_ = 0.108;
        beta_ = 0.3744;
        R_ = d_max_ / d_min_;
        Fai_ = 0.6;
        C_ = 0.2712;
        epsilon_ = 2.0957;
        gravity_ = 9.8;
      }
//=============================================================================================//
void UpdateDiffusivityAndSegregationRate::update(size_t index_i, Real dt)
{
    Real concentration_i = concentration_[index_i];
    Real d_avg = concentration_i * d_min_ + (1-concentration_i) * d_max_;
    test_[index_i] = d_avg;
    //Real d_avg = 0.3;
    Real gama_i = gama_[index_i];
    Real rho_i = rho_[index_i];
    Real p_i = SMIN(p_[index_i], 1000*10*0.01);
    //Real p_i = p_[index_i];
    /*Miu(I) model parameters*/
    // inertial_num_[index_i] = d_avg* gama_i/sqrt(p_i/rho_[index_i] + TinyReal);
    // Real miu_i = plastic_continuum_.getViscosityMiuI(inertial_num_[index_i], 0.6,1.0);
    // alpha_each_[index_i] = plastic_continuum_.getDPConstantsA_WithMiu(miu_i);
    // Kc_each_[index_i] = plastic_continuum_.getDPConstantsK_WithMiu(miu_i);
    
    //---------------------------/*Gray et al 2021*/-------------------------//
    // /*Segregation Velocity*/
    // Real F_ = R_ - 1 + epsilon_ * (1-concentration_i)*(R_-1)*(R_-1);
    // Real upper_term = beta_  * rho_i * gravity_* gama_i * d_avg * d_avg * F_ ;
    // Real lower_term =  C_ * rho_i * gravity_* d_avg + p_i;
    // Real segreation_rate_i = upper_term / lower_term;
    // segregation_rate_[index_i] = segreation_rate_i;
    // /*Diffusion Rate*/
    // Real diffusion_rate_i = A_ * gama_i * d_avg * d_avg ;
    // diffusion_rate_[index_i] = diffusion_rate_i;

    //---------------------------/*Zhu et al*/-------------------------//
    segregation_rate_[index_i] = x_v_ * gama_i * d_avg;
    diffusion_rate_[index_i] = x_d_ * gama_i * d_avg * d_avg;
}
//=============================================================================================//
UpdateDiffusivityAndSegregationRateGray::UpdateDiffusivityAndSegregationRateGray(SPHBody &sph_body, Real height, Real period)
    : LocalDynamics(sph_body),
      concentration_(particles_->getVariableDataByName<Real>("Concentration")),
      gama_(particles_->getVariableDataByName<Real>("ShearRate")),
      segregation_rate_(particles_->getVariableDataByName<Real>("SegregationRate")),
      diffusion_rate_(particles_->getVariableDataByName<Real>("SegregationDiffusivity")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      height_(height), period_(period)
      {
        d_min_ = 0.004;
        d_max_ = 0.008;
        A_ = 0.108;
        beta_ = 0.3744;
        R_ = d_max_ / d_min_;
        Fai_ = 0.6;
        C_ = 0.2712;
        epsilon_ = 2.0957;
        actual_height_ = 1.0;
      }
//=============================================================================================//
void UpdateDiffusivityAndSegregationRateGray::update(size_t index_i, Real dt)
{
    Real concentration_i = concentration_[index_i];
    //Real concentration_i =0.4;
    Real d_avg = concentration_i * d_min_ + (1-concentration_i) * d_max_;
    //Real d_avg = 0.006;
    // Real gama_i = gama_[index_i];
    Real gama_i = gama_[index_i] * 10 / 0.3;
    /*Segregation Velocity*/
    Real F_ = R_ - 1 + epsilon_ * (1-concentration_i)*(R_-1)*(R_-1);
    //Real F_ = R_ - 1;
    Real pos_y = pos_[index_i][1];
    Real gama_t = 2.308;
    //Real gama_t = gama_i * period_;
    Real upper_term = beta_  * gama_t * d_avg * d_avg * F_ / height_ /height_;
    //Real lower_term = C_ * d_avg / height_ + Fai_ * (1 - pos_y/actual_height_);
    Real relative_height = pos_y/actual_height_;
    relative_height = SMIN(0.9,relative_height);
    //relative_height = SMAX(0.1,relative_height);
    Real lower_term = Fai_ * (1 - relative_height) +  C_ * d_avg / height_;
    Real segreation_rate_i = upper_term / lower_term;
    segregation_rate_[index_i] = upper_term /lower_term ;

    /*Diffusion Rate*/
    Real diffusion_rate_i = A_ * gama_t * d_avg * d_avg / height_ / height_;
    diffusion_rate_[index_i] = diffusion_rate_i;
    
}
//=============================================================================================//
UpdateDiffusivityAndSegregationRateGray2023NonDim::UpdateDiffusivityAndSegregationRateGray2023NonDim(SPHBody &sph_body, Real height, Real period)
    : LocalDynamics(sph_body),
      concentration_(particles_->getVariableDataByName<Real>("Concentration")),
      gama_(particles_->getVariableDataByName<Real>("ShearRate")),
      segregation_rate_(particles_->getVariableDataByName<Real>("SegregationRate")),
      diffusion_rate_(particles_->getVariableDataByName<Real>("SegregationDiffusivity")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      height_(height), period_(period)
      {
        d_mean_ = 0.006;
        d_min_ = 0.004;
        d_max_ = 0.008;
        A_ = 0.108;
        beta_ = 0.3744;
        R_ = d_max_ / d_min_;
        Fai_ = 0.6;
        C_ = 0.2712;
        epsilon_ = 2.0957;
        actual_height_ = 0.09 / d_mean_;
        U_ = sqrt(10.0*0.006);
      }
//=============================================================================================//
void UpdateDiffusivityAndSegregationRateGray2023NonDim::update(size_t index_i, Real dt)
{
    Real concentration_i = concentration_[index_i];
    //Real concentration_i =0.4;
    Real d_avg = concentration_i * d_min_ + (1-concentration_i) * d_max_;
    d_avg = d_avg/d_mean_;
    //Real d_avg = 0.006;
    Real gama_i = gama_[index_i];
    /*Segregation Velocity*/
    Real F_ = R_ - 1 + epsilon_ * (1-concentration_i)*(R_-1)*(R_-1);
    //Real F_ = R_ - 1;
    Real pos_y = pos_[index_i][1]/d_mean_;
    //Real gama_t = 2.308 * d_mean_ / U_;
    // Real T_ = 13.0 * U_ / d_mean_;
    // gama_i = 2.308 / T_;
    Real upper_term = beta_  * gama_i * d_avg * d_avg * F_;
    //Real lower_term = C_ * d_avg / height_ + Fai_ * (1 - pos_y/actual_height_);
    Real relative_height = actual_height_ - pos_y ;
    // relative_height = SMIN(0.9 * actual_height_,relative_height);
    // relative_height = SMAX(0.1 * actual_height_,relative_height);
    //Real relative_height = 0.5 * actual_height_;
    Real lower_term = Fai_ * relative_height+  C_ * d_avg;
    Real segreation_rate_i = upper_term / lower_term;
    segregation_rate_[index_i] = segreation_rate_i;

    /*Diffusion Rate*/
    Real diffusion_rate_i = A_ * gama_i * d_avg * d_avg ;
    diffusion_rate_[index_i] = diffusion_rate_i;
    
}
//=============================================================================================//
UpdateDiffusivityAndSegregationRateGray2021::UpdateDiffusivityAndSegregationRateGray2021(SPHBody &sph_body, Real height, Real period)
    : LocalDynamics(sph_body),
      concentration_(particles_->getVariableDataByName<Real>("Concentration")),
      gama_(particles_->getVariableDataByName<Real>("ShearRate")),
      segregation_rate_(particles_->getVariableDataByName<Real>("SegregationRate")),
      diffusion_rate_(particles_->getVariableDataByName<Real>("SegregationDiffusivity")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      height_(height), period_(period)
      {
        d_min_ = 0.5;
        d_max_ = 0.1;
        A_ = 0.108;
        beta_ = 0.3744;
        R_ = d_max_ / d_min_;
        Fai_ = 0.6;
        C_ = 0.2712;
        epsilon_ = 2.0957;
        gravity_ = 9.8;
      }
//=============================================================================================//
void UpdateDiffusivityAndSegregationRateGray2021::update(size_t index_i, Real dt)
{
    Real concentration_i = concentration_[index_i];
    Real d_avg = concentration_i * d_min_ + (1-concentration_i) * d_max_;
    Real gama_i = gama_[index_i];
    Real rho_i = rho_[index_i];
    Real pos_y = pos_[index_i][1];
    Real p_i = p_[index_i];
    //Real p_i = 0.6 * rho_i * gravity_ * SMAX(height_-pos_y, 0.3*height_);
    //Real p_i = rho_i * gravity_ * height_ * 0.5;
    /*Segregation Velocity*/
    Real F_ = R_ - 1 + epsilon_ * (1-concentration_i)*(R_-1)*(R_-1);
    Real upper_term = beta_  * rho_i * gravity_* gama_i * d_avg * d_avg * F_ ;
    Real lower_term =  C_ * rho_i * gravity_* d_avg + p_i;
    Real segreation_rate_i = upper_term / lower_term;
    //segreation_rate_i *= period_ / height_;
    segregation_rate_[index_i] = segreation_rate_i;

    /*Diffusion Rate*/
    Real diffusion_rate_i = A_ * gama_i * d_avg * d_avg ;
    //diffusion_rate_i *= period_ / height_ / height_ 
    diffusion_rate_[index_i] = diffusion_rate_i;
    
}
//=================================================================================================//
} // namespace SPH
