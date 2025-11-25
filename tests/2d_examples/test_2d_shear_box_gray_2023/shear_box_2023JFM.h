/**
 * @file 	segregation_ck.h
 * @brief 	This is the case file for the test of particle segregation
 * @author  Shaung Li and Xiangyu Hu
 */
#ifndef SHEAR_BOX_JFM2023_H
#define SHEAR_BOX_JFM2023_H
#include "sphinxsys.h"
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real mean_d = 0.006;
Real U_ = sqrt(0.006*10);
Real LL = 0.04/mean_d;                       /**< Soil column length. */
Real LH = 0.09/mean_d;                       /**< Soil column height. */
Real DL = LL;                       /**< Tank length. */
Real DH = 1.0 * LH;                      /**< Tank height. */
Real particle_spacing_ref = LH/50; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;  /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
//----------------------------------------------------------------------
//	Material properties of the soil.
//----------------------------------------------------------------------
Real rho0_s = 2040;                                                       // reference density of soil
Real gravity_g = 9.8;                                                     // gravity force of soil
Real Youngs_modulus = 5.84e6;                                             // reference Youngs modulus
Real poisson = 0.3;                                                       // Poisson ratio
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3.0 * (1.0 - 2.0 * poisson))); // sound speed
Real friction_angle = 21.9 * Pi / 180;
Real cohesion = 0.0;
Real dilatancy = 0.0;
 
//Diffusion
Real shear_period = 1.0  * U_/ mean_d;
std::string diffusion_species_name = "Concentration";
Real B = LH;
Real Sr = 0.03;
Real Dr = 0.03/29.6;
Real v = Sr * 0.09 /1.0;
Real d = Dr * 0.09 * 0.09 /1.0;
Real segregation_rate = v/U_;
Real diffusion_coeff =  d/U_/mean_d;
Real initial_temperature = 0.5;
Real heat_flux = 0.0; //from the Nemann boundary
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d soil_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin:
Vec2d soil_block_translation = soil_block_halfsize;
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
Vec2d lifted_outer_wall_translation = inner_wall_halfsize + Vecd(0.0,BW); 
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


/*Shear part*/
MultiPolygon createSoilShearShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addABox(Transform(soil_block_translation), soil_block_halfsize, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createShearShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addABox(Transform(lifted_outer_wall_translation), outer_wall_halfsize, ShapeBooleanOps::add);

    return multi_polygon;
}
//----------------------------------------------------------------------//
/*Shear dynamics*/
//----------------------------------------------------------------------//
class ShearMaking : public BodyPartMotionConstraint
{
    Real gama_0_;
    Real period_;
    Real omega_;

    Vecd getDisplacement(const Real &time, const Real &height)
    {
        Vecd displacement{Vecd::Zero()};
        displacement[0] = height * gama_0_ * sin(omega_ * time);
        return displacement;
    }

    Vec2d getVelocity(const Real &time, const Real &height)
    {
        Vec2d velocity{Vecd::Zero()};
        velocity[0] = height * gama_0_ * omega_ * cos(omega_ * time);
        return velocity;
    }

    Vec2d getAcceleration(const Real &time, const Real &height)
    {
        Vec2d acceleration{Vecd::Zero()};
        acceleration[0] = -1.0 * gama_0_ * height * omega_ * omega_ * sin(omega_ * time);
        return acceleration;
    }


  public:
    ShearMaking(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part),
          gama_0_(0.577),
          period_(shear_period),
          acc_(particles_->registerStateVariable<Vecd>("Acceleration")),
          physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime"))
    {
      omega_ = 2.0 * PI / period_;
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        Real time = * physical_time_;
        Real height = pos0_[index_i][1];
        pos_[index_i] = pos0_[index_i] + getDisplacement(time, height);
        vel_[index_i] = getVelocity(time, height);
        acc_[index_i] = getAcceleration(time, height);
    };

  protected:
    Vecd *acc_;
    Real *physical_time_;
};


/*Test new shear - consider vertical displacement*/
// class ShearMaking : public BodyPartMotionConstraint
// {
//     Real gama_0_;
//     Real period_;
//     Real omega_;

//     Vecd getDisplacement(const Real &time, const Real &height)
//     {
//         Vecd displacement{Vecd::Zero()};
//         displacement[0] = height * gama_0_ * sin(omega_ * time);
//         displacement[1] = -0.134 * height * ABS(sin(omega_ * time)); // 添加体积守恒垂向位移
//         //displacement[1] = 0.0; // 添加体积守恒垂向位移
//         return displacement;
//     }

//     Vec2d getVelocity(const Real &time, const Real &height)
//     {
//         Vec2d velocity{Vecd::Zero()};
//         velocity[0] = height * gama_0_ * omega_ * cos(omega_ * time);
//         //velocity[1] = -0.134 * height * omega_ * cos(omega_ * time); // 添加体积守恒垂向速度
//         return velocity;
//     }

//     Vec2d getAcceleration(const Real &time, const Real &height)
//     {
//         Vec2d acceleration{Vecd::Zero()};
//         acceleration[0] = -1.0 * gama_0_ * height * omega_ * omega_ * sin(omega_ * time);
//         //acceleration[1] = -0.134 * height * omega_ * omega_* sin(omega_ * time); // 添加体积守恒垂向加速度
//         return acceleration;
//     }

// public:
//     ShearMaking(BodyPartByParticle &body_part)
//         : BodyPartMotionConstraint(body_part),
//           gama_0_(0.577),
//           period_(shear_period),
//           acc_(particles_->registerStateVariable<Vecd>("Acceleration")),
//           physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime"))
//     {
//         omega_ = 2.0 * PI / period_;
//     }

//     void update(size_t index_i, Real dt = 0.0)
//     {
//         Real time = *physical_time_;
//         Real height = pos0_[index_i][1];
//         pos_[index_i] = pos0_[index_i] + getDisplacement(time, height);
//         vel_[index_i] = getVelocity(time, height);
//         acc_[index_i] = getAcceleration(time, height);
//     }

// protected:
//     Vecd *acc_;
//     Real *physical_time_;
// };


class SheaboxDiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit SheaboxDiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)),
          pos_(particles_->registerStateVariable<Vecd>("Position")){};

    void update(size_t index_i, Real dt)
    {
        if(pos_[index_i][1]>0.5 * LH)
            phi_[index_i] = 1.0;
        else
            phi_[index_i] = 0.0;
        // phi_[index_i] = pos_[index_i][1];
    };

  protected:
    Real *phi_;
    Vecd *pos_;
};

class LinearSheaboxDiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit LinearSheaboxDiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)),
          pos_(particles_->registerStateVariable<Vecd>("Position")){};

    void update(size_t index_i, Real dt)
    {
      phi_[index_i] = pos_[index_i][1];
    };

  protected:
    Real *phi_;
    Vecd *pos_;
};

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
#endif