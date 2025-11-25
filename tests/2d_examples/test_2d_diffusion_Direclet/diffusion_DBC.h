#pragma once
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 0.02;
Real H = 0.002;
Real resolution_ref = H / 25.0;
Real BW = resolution_ref * 4.0;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1.43e-7;
std::string diffusion_species_name = "Phi";
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 20.0;
Real phi_upper_wall = 20.0;
Real phi_lower_wall = 40.0;
Real phi_fluid_initial = 20.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
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

std::vector<Vecd> lower_region_edge_points{
    Vecd(-BW, -BW), Vecd(-BW, 0.0), Vecd(L+BW, 0.0),
    Vecd(L+BW, -BW), Vecd(-BW, -BW)};
std::vector<Vecd> upper_region_edge_points{
    Vecd(-BW, H), Vecd(-BW, H + BW), Vecd(L+BW, H + BW),
    Vecd(L+BW, H), Vecd(-BW, H)};
//----------------------------------------------------------------------
// Define extra classes which are used in the main program.
// These classes are defined under the namespace of SPH.
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
class DiffusionBody : public MultiPolygonShape
{
  public:
    explicit DiffusionBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::add);
    }
};

class DirichletWallBoundary : public MultiPolygonShape
{
  public:
    explicit DirichletWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(lower_region_edge_points, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(upper_region_edge_points, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
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

//----------------------------------------------------------------------
//	Specify diffusion relaxation method.
//----------------------------------------------------------------------
using DiffusionBodyRelaxation = DiffusionBodyRelaxationComplex<
    IsotropicDiffusion, KernelGradientInner, KernelGradientContact, Dirichlet>;

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;

    observation_points.push_back(Vecd(0.5 * L, 0.5 * H));

    return observation_points;
};
