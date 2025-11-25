/**
 * @file 	segregation_ck.h
 * @brief 	This is the case file for the test of particle segregation
 * @author  Shaung Li and Xiangyu Hu
 */
#ifndef SHEAR_BOX_H
#define SHEAR_BOX_H
#include "sphinxsys.h"
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 0.02;
Real H = 0.002;
Real particle_spacing_ref = H/25.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;  /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
// observer location
StdVec<Vecd> observation_location = {Vecd(0.5 * L, 0.5 * H)};
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
Real shear_period = 1.0;
std::string diffusion_species_name = "Concentration";
Real segregation_rate = 1.43e-7;
Real diffusion_coeff = 1.43e-7;
//----------------------------------------------------------------------
Real initial_temperature = 20.0;
Real phi_upper_wall = 20.0;
Real phi_lower_wall = 40.0;
Real phi_fluid_initial = 20.0;
Real heat_flux = 0.0; //from the Nemann boundary

//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
std::vector<Vecd> lower_region_edge_points{
    Vecd(-BW, -BW), Vecd(-BW, 0.0), Vecd(L+BW, 0.0),
    Vecd(L+BW, -BW), Vecd(-BW, -BW)};
std::vector<Vecd> upper_region_edge_points{
    Vecd(-BW, H), Vecd(-BW, H + BW), Vecd(L+BW, H + BW),
    Vecd(L+BW, H), Vecd(-BW, H)};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(lower_region_edge_points, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(upper_region_edge_points, ShapeBooleanOps::add);
    }
};


std::vector<Vecd> createThermalDomain()
{
    std::vector<Vecd> thermalDomainShape;
    thermalDomainShape.push_back(Vecd(0.0, 0.0));
    thermalDomainShape.push_back(Vecd(0.0, H));
    thermalDomainShape.push_back(Vecd(L, H));
    thermalDomainShape.push_back(Vecd(L, 0.0));
    thermalDomainShape.push_back(Vecd(0.0, 0.0));

    return thermalDomainShape;
}

class SoilBlock : public MultiPolygonShape
{
  public:
    explicit SoilBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------//
class DiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit DiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature;
    };

  protected:
    Real *phi_;
};

class DirichletWallBoundaryInitialCondition : public LocalDynamics
{
  public:
    explicit DirichletWallBoundaryInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariable<Real>(diffusion_species_name)){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = -0.0;

        if (pos_[index_i][1] > 0.5 * H)
        {
            phi_[index_i] = phi_upper_wall;
        }
        if (pos_[index_i][1] < 0.5 *H )
        {
            phi_[index_i] = phi_lower_wall;
        }
    }

  protected:
    Vecd *pos_;
    Real *phi_;
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;

    observation_points.push_back(Vecd(0.5 * L, 0.5 * H));

    return observation_points;
};


//----------------------------------------------------------------------//
using DiffusionBodyRelaxation = SegregationDiffusionBodyRelaxationComplex<
    SegregationLocalIsotropicDiffusion, KernelGradientInner, KernelGradientContact, Dirichlet>;

using CorrectedDiffusionBodyRelaxation = SegregationDiffusionBodyRelaxationComplex<
    SegregationLocalIsotropicDiffusion, CorrectedKernelGradientInner, CorrectedKernelGradientContact, Dirichlet>;
#endif