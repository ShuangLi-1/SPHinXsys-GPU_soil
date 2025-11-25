#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;
#define PI 3.1415926
// general parameters for geometry
Real radius = 0.1;                                         // Soil column length
Real height = 0.1;                                         // Soil column height
Real resolution_ref = 0.01;                         // particle spacing
Real BW = resolution_ref * 4;                              // boundary width
Real DL = 3.0; // tank length
Real DH = 3.0;                                   // tank height
Real DW = 3.0;                                          // tank width
// for material properties
Real rho0_s = 2040;                                                       // reference density of soil
Real gravity_g = 9.8;                                                     // gravity force of soil
Real Youngs_modulus = 5.84e6;                                             // reference Youngs modulus
Real poisson = 0.3;                                                       // Poisson ratio
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3 * (1 - 2 * poisson)));
Real friction_angle = 30.0 * Pi / 180;
Real cohesion = 3000.0;
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

std::string soil_full_path_to_file = "./input/SlideBody.stl";  //4+45 有incline
std::string wall_side_full_path_to_file = "./input/Channel.stl";

class SoilBlock : public ComplexShape
{
  public:
    explicit SoilBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(soil_full_path_to_file, Vecd(0.00, 0.0, 0.0), 1.0);
    }
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(wall_side_full_path_to_file, Vecd(0, 0.0, 0), 1.0);
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