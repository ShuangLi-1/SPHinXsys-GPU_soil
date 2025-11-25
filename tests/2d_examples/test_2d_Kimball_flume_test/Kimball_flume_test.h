/**
 * @file 	segregation_ck.h
 * @brief 	This is the case file for the test of particle segregation
 * @author  Shaung Li and Xiangyu Hu
 */
#pragma once
#include "sphinxsys.h"
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 20.0;                       /**< Tank length. */
Real DH = 5;                      /**< Tank height. */
Real LL = 0.5;                       /**< Soil column length. */
Real LH = 0.3;                       /**< Soil column height. */
Real particle_spacing_ref = 0.025; /**< Initial reference particle spacing. */
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
Real c_s = 30; // sound speed
Real friction_angle = 90.0 * Pi / 180;
Real cohesion = 2000.0;
Real dilatancy = 0.0;
 
//Diffusion
Real shear_period = 1.0;
std::string diffusion_species_name = "Concentration";
Real B = 1.0;
Real Sr = 0.016 * 20.0;
Real Dr = Sr/20.9;
Real segregation_rate = Sr * B / shear_period;
Real diffusion_coeff =  Dr * B * B / shear_period;
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

Real x_offset = 5.0;
Real y_offset = -x_offset * tan(30.0 * Pi / 180);
//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallBlockShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(-BW, -BW));
    pnts.push_back(Vecd(-BW, 4.115));
    pnts.push_back(Vecd(7.13, 0.0));
    pnts.push_back(Vecd(17.0, 0.0));
    pnts.push_back(Vecd(17.0, -BW));
    pnts.push_back(Vecd(-BW, -BW));

    return pnts;
}

std::vector<Vecd> createInnerWallBlockShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(-BW, -BW));
    pnts.push_back(Vecd(-BW, 4.115- BW));
    pnts.push_back(Vecd(7.13 - BW, -BW));
    pnts.push_back(Vecd(-BW, -BW));
    return pnts;
}

std::vector<Vecd> createWallBlockShapeForSub()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(-BW, -BW));
    pnts.push_back(Vecd(-BW, 4.115));
    pnts.push_back(Vecd(7.13, 0.0));
    pnts.push_back(Vecd(-BW, -BW));
    return pnts;
}
std::vector<Vecd> createRightWallShape()
{
    std::vector<Vecd> pnts;
    Real wall_pos = 9.0;
    pnts.push_back(Vecd(wall_pos, -BW));
    pnts.push_back(Vecd(wall_pos, 2.0));
    pnts.push_back(Vecd(wall_pos+BW, 2.0));
    pnts.push_back(Vecd(wall_pos+BW, -BW));
    pnts.push_back(Vecd(wall_pos, -BW));
    return pnts;
}

/*original*/
// std::vector<Vecd> createSoilBlockShape()
// {
//     std::vector<Vecd> pnts;
//     pnts.push_back(Vecd(0.0, 3.37));
//     pnts.push_back(Vecd(0.0, 4.115));
//     pnts.push_back(Vecd(1.3, 4.115));
//     pnts.push_back(Vecd(1.3, 3.37)); // lower right
//     pnts.push_back(Vecd(0.0, 3.37));
//     return pnts;
// }

/*Test*/
std::vector<Vecd> createSoilBlockShape()
{
    std::vector<Vecd> pnts;
    Vecd offset(x_offset, y_offset);
    pnts.push_back(Vecd(0.0, 3.37)+ offset) ;
    pnts.push_back(Vecd(0.0, 4.115)+ offset) ;
    pnts.push_back(Vecd(1.3, 4.115)+ offset) ;
    pnts.push_back(Vecd(1.3, 3.37)+ offset) ; // lower right
    pnts.push_back(Vecd(0.0, 3.37)+ offset) ;
    return pnts;
}


class WallBoundary : public MultiPolygonShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallBlockShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createRightWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallBlockShape(), ShapeBooleanOps::sub);
    }
};



class SoilBlock : public MultiPolygonShape
{
public:
    explicit SoilBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createSoilBlockShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createWallBlockShapeForSub(), ShapeBooleanOps::sub);
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
    Vec2d lifted_outer_wall_translation = inner_wall_halfsize + Vecd(0.0,BW); 
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
        Real time = *physical_time_;
        Real height = pos0_[index_i][1];
        pos_[index_i] = pos0_[index_i] + getDisplacement(time, height);
        vel_[index_i] = getVelocity(time, height);
        acc_[index_i] = getAcceleration(time, height);
    };

  protected:
    Vecd *acc_;
    Real *physical_time_;
};

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