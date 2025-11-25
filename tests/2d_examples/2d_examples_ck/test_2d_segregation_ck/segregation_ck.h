/**
 * @file 	segregation_ck.h
 * @brief 	This is the case file for the test of particle segregation
 * @author  Shaung Li and Xiangyu Hu
 */
#ifndef TEST_2D_SEGREGATION_H
#define TEST_2D_SEGREGATION_H
#include "sphinxsys_ck.h"
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.04;                       /**< Tank length. */
Real DH = 0.15;                      /**< Tank height. */
Real LL = 0.04;                       /**< Soil column length. */
Real LH = 0.09;                       /**< Soil column height. */
Real particle_spacing_ref = LH / 25; /**< Initial reference particle spacing. */
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
std::string diffusion_species_name = "Concentration";
Real diffusion_coeff = 2.07e-7;
Real segregation_rate = 1.2e-4;
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

class SoilBlock : public MultiPolygonShape
{
  public:
    explicit SoilBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(soil_shape, ShapeBooleanOps::add);
    }
};


/*Shear part*/
MultiPolygon createShearShape()
{
    std::vector<Vecd> wave_make_shape;
    wave_make_shape.push_back(Vecd(-BW, 0.0));
    wave_make_shape.push_back(Vecd(-BW, DH + BW));
    wave_make_shape.push_back(Vecd(0.0, DH + BW));
    wave_make_shape.push_back(Vecd(0.0, 0.0));
    wave_make_shape.push_back(Vecd(-BW, 0.0));

    MultiPolygon multi_polygon;
    //multi_polygon.addAPolygon(wave_make_shape, ShapeBooleanOps::add);
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
          period_(13),
          acc_(particles_->registerStateVariable<Vecd>("Acceleration")),
          physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime"))
    {
      omega_ = 2.0 * PI / period_;
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        Real time = *physical_time_;
        Real height = pos_[index_i][1];
        pos_[index_i] = pos0_[index_i] + getDisplacement(time, height);
        vel_[index_i] = getVelocity(time, height);
        acc_[index_i] = getAcceleration(time, height);
    };

  protected:
    Vecd *acc_;
    Real *physical_time_;
};

Real h = 1.3 * particle_spacing_ref;
MultiPolygon createWaveProbeShape4()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(3.99 - h, 0.0));
    pnts.push_back(Vecd(3.99 - h, 1.0));
    pnts.push_back(Vecd(3.99 + h, 1.0));
    pnts.push_back(Vecd(3.99 + h, 0.0));
    pnts.push_back(Vecd(3.99 - h, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}
MultiPolygon createWaveProbeShape5()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(7.02 - h, 0.155));
    pnts.push_back(Vecd(7.02 - h, 1.0));
    pnts.push_back(Vecd(7.02 + h, 1.0));
    pnts.push_back(Vecd(7.02 + h, 0.155));
    pnts.push_back(Vecd(7.02 - h, 0.155));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}

MultiPolygon createWaveProbeShape12()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(8.82 - h, 0.155));
    pnts.push_back(Vecd(8.82 - h, 1.0));
    pnts.push_back(Vecd(8.82 + h, 1.0));
    pnts.push_back(Vecd(8.82 + h, 0.155));
    pnts.push_back(Vecd(8.82 - h, 0.155));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}

StdVec<Vecd> creatObserverPositions()
{
    StdVec<Vecd> observer_positions;
    observer_positions.push_back(Vecd(7.862, 0.645));
    observer_positions.push_back(Vecd(7.862, 0.741));
    observer_positions.push_back(Vecd(7.862, 0.391));
    observer_positions.push_back(Vecd(7.862, 0.574));
    observer_positions.push_back(Vecd(7.862, 0.716));
    observer_positions.push_back(Vecd(7.862, 0.452));
    return observer_positions;
}
//----------------------------------------------------------------------//
/*Thermal conditions and functions*/
//----------------------------------------------------------------------//
class ThermosolidBodyInitialCondition : public LocalDynamics
{
  public:
    explicit ThermosolidBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {

        phi_[index_i] = 0.5;

    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------//
using ThermalRelaxationComplex = DiffusionBodyRelaxationComplex<
    IsotropicDiffusion, KernelGradientInner, KernelGradientContact, Dirichlet>;

#endif