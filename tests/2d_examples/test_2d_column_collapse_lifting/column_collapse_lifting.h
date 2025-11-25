/**
 * @file 	column_collapse.cpp
 * @brief 	2D column collapse.
 * @details This is the one of the basic test cases, also the first case for understanding
 * 			SPH method for modelling granular materials such as soils and sands.
 * @author Shuaihao Zhang and Xiangyu Hu
 */
#pragma once
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.4;                       /**< Tank length. */
Real DH = 0.15;                      /**< Tank height. */
Real LL = 0.08;                       /**< Soil column length. */
Real LH = 0.05;                       /**< Soil column height. */
Real GL = LL;                       /**< Gate length. */
Real GH = LH;                       /**< Gate height. */
Real particle_spacing_ref = LH / 60; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;  /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the soil.
//----------------------------------------------------------------------
Real rho0_s = 2040;                                                       // reference density of soil
Real gravity_g = 9.8;                                                     // gravity force of soil
Real Youngs_modulus = 5.84e6;                                             // reference Youngs modulus
Real poisson = 0.3;                                                       // Poisson ratio
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3.0 * (1.0 - 2.0 * poisson))); // sound speed
Real friction_angle = 32.5 * Pi / 180;
Real cohesion = 0.0;
Real dilatancy = 0.0;
Real miu_s = tan(20.9 * Pi / 180);
Real miu_d = tan(32.5 * Pi / 180);
Real d_min = 0.002;
Real d_max = 0.002;
//Diffusion
Real shear_period = 1.0;
std::string diffusion_species_name = "Concentration";
Real B = 1.0;
Real Sr = 0.016 * 5.0;
Real Dr = Sr/20.9;
Real segregation_rate = Sr * B / shear_period;
Real diffusion_coeff =  100.0 * Dr * B * B / shear_period;
Real initial_temperature = 0.5;
Real heat_flux = 0.0; //from the Nemann boundary
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d soil_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin:
Vec2d soil_block_translation = soil_block_halfsize + Vecd(0.5 * DL,0.0); 
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;

Vec2d outer_gate_halfsize = Vecd(0.5 * LL + BW, 0.5 * LH);
Vec2d outer_gate_translation = soil_block_translation;
Vec2d inner_gate_halfsize = Vecd(0.5 * LL, 0.5 * LH);
Vec2d inner_gate_translation = soil_block_translation;

//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
        add<TransformShape<GeometricShapeBox>>(Transform(outer_gate_translation), outer_gate_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_gate_translation), inner_gate_halfsize);
    }
};
std::vector<Vecd> soil_shape{
    Vecd(0, 0), Vecd(0, LH), Vecd(LL, LH), Vecd(LL, 0), Vecd(0, 0)};

class SoilBlock : public ComplexShape
{
  public:
    explicit SoilBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
      add<TransformShape<GeometricShapeBox>>(Transform(soil_block_translation), soil_block_halfsize);
    }
};

//----------------------------------------------------------------------//
class SheaboxDiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit SheaboxDiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)),
          pos_(particles_->registerStateVariable<Vecd>("Position")){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = 0.5;
    };


  protected:
    Real *phi_;
    Vecd *pos_;
};
//----------------------------------------------------------------------//

class NeumannWallBoundaryInitialCondition : public LocalDynamics
{
  public:
    explicit NeumannWallBoundaryInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)),
          phi_flux_(particles_->getVariableDataByName<Real>(diffusion_species_name + "Flux")) {}

    void update(size_t index_i, Real dt)
    {
        phi_flux_[index_i] = 0.0;
    }

  protected:
    Vecd *pos_;
    Real *phi_, *phi_flux_;
};
//----------------------------------------------------------------------//
using DiffusionBodyRelaxation = SegregationDiffusionBodyRelaxationComplex<
    SegregationLocalIsotropicDiffusion, KernelGradientInner, KernelGradientContact, Neumann>;

using CorrectedDiffusionBodyRelaxation = SegregationDiffusionBodyRelaxationComplex<
    SegregationLocalIsotropicDiffusion, CorrectedKernelGradientInner, CorrectedKernelGradientContact, Neumann>;

/*Lifting shape*/
MultiPolygon createLiftingShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addABox(Transform(outer_gate_translation), outer_gate_halfsize, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------//
/*Lifting dynamics*/
//----------------------------------------------------------------------//
class Lifting : public BodyPartMotionConstraint
{
    Real a_up_; // 向上加速度，单位 m/s^2

    Vecd getDisplacement(const Real &time) const
    {
        Vecd displacement = Vecd::Zero();
        displacement[1] = SMIN(2.0*DH, 0.5 * a_up_ * time * time); // y方向位移
        return displacement;
    }

    Vecd getVelocity(const Real &time) const
    {
        Vecd velocity = Vecd::Zero();
        velocity[1] = a_up_ * time; // y方向速度
        return velocity;
    }

    Vecd getAcceleration() const
    {
        Vecd acceleration = Vecd::Zero();
        acceleration[1] = a_up_; // y方向恒定加速度
        return acceleration;
    }

  public:
    explicit Lifting(BodyPartByParticle &body_part, Real a_up = 8.0)
        : BodyPartMotionConstraint(body_part),
          a_up_(a_up),
          acc_(particles_->registerStateVariable<Vecd>("Acceleration")),
          physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime"))
    {}

    void update(size_t index_i, Real dt = 0.0) 
    {
        Real time = *physical_time_;
        pos_[index_i] = pos0_[index_i] + getDisplacement(time);
        vel_[index_i] = getVelocity(time);
        acc_[index_i] = getAcceleration();
    }

  protected:
    Vecd *acc_;
    Real *physical_time_;
};

