#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;
#define PI 3.1415926
// general parameters for geometry
Real radius = 100;                                         // Soil column length
Real height = 100;                                         // Soil column height
Real resolution_ref = 10.0;                         // particle spacing
Real BW = resolution_ref * 4;                              // boundary width
Real DL = radius; // tank length
Real DH = 1000;                                   // tank height
Real DW = DL;                                              // tank width
// for material properties
Real rho0_s = 1540;           // reference density of soil
Real gravity_g = 9.8;         // gravity force of soil
Real Youngs_modulus = 10.0e6; // reference Youngs modulus
Real poisson = 0.2;           // Poisson ratio
Real U_max_s = 80; //By test
//Real c_s = sqrt(Youngs_modulus / (rho0_s * 3 * (1 - 2 * poisson)));
Real c_s = 500;
Real friction_angle = 10.5* Pi / 180;
Real cohesion = 15000;
Real dilatancy = 0.0;

//Real U_max_s  = 2.0 * sqrt(gravity_g * 0.15);; //By test
/** Define the soil body. */
Real inner_circle_radius = radius;
int resolution(20/100);

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

std::string soil_full_path_to_file = "./input/chain-soil.stl";  //4+45 有incline
std::string wall_side_full_path_to_file = "./input/chain-wall-side.stl";
std::string wall_base_full_path_to_file = "./input/chain-wall-base.stl";

class SoilBlock : public ComplexShape
{
  public:
    explicit SoilBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column(DL / 2, 0.5 * height, DW / 2);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), inner_circle_radius,
            0.5 * height, resolution, translation_column);

    }
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
      Vecd translation_column(DL / 2, 0.5 * height, DW / 2);
      add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), inner_circle_radius,
            0.5 * height, resolution, translation_column);


    }
};


class DiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit DiffusionInitialCondition(SPHBody &sph_body)
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