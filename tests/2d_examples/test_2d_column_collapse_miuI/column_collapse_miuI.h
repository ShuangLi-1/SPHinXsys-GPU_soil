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
Real DL = 0.5;                       /**< Tank length. */
Real DH = 0.15;                      /**< Tank height. */
Real LL = 0.2;                       /**< Soil column length. */
Real LH = 0.1;                       /**< Soil column height. */
Real particle_spacing_ref = LH / 25; /**< Initial reference particle spacing. */
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
Real friction_angle = 21.9 * Pi / 180;
Real cohesion = 0.0;
Real dilatancy = 0.0;
//Diffusion
Real shear_period = 1.0;
std::string diffusion_species_name = "Concentration";
Real B = 1.0;
Real Sr = 0.016 * 5.0;
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
