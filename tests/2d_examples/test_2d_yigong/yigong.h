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
Real LL = 10000;                       /**< Soil column length. */
Real LH = 4000;                       /**< Soil column height. */
Real DL = LL;                       /**< Tank length. */
Real DH = 1.0 * LH;                      /**< Tank height. */
Real particle_spacing_ref = 10.0; /**< Initial reference particle spacing. */
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
Real c_s = 500; // sound speed
Real friction_angle = 30.0 * Pi / 180;
Real cohesion = 30e3;
Real dilatancy = 0.0;
 
//Diffusion
Real shear_period = 1.0;
std::string diffusion_species_name = "Concentration";
Real B = 1.0;
Real Sr = 0.016 * 20;
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
std::vector<Vecd> createWallShape()
{
  std::vector<Vecd> pnts;
  pnts.push_back(Vecd(32.96250515039128, 3213.174748398902));
  pnts.push_back(Vecd(98.8875154511743, 3139.981701738335));
  pnts.push_back(Vecd(370.8281829419034, 2766.6971637694414));
  pnts.push_back(Vecd(824.0626287597854, 2268.9844464775842));
  pnts.push_back(Vecd(1466.8314791924186, 1910.3385178408053));
  pnts.push_back(Vecd(2216.7284713638232, 1602.9277218664229));
  pnts.push_back(Vecd(2884.2192006592495, 1302.8362305580968));
  pnts.push_back(Vecd(3683.559950556242, 1090.5763952424513));
  pnts.push_back(Vecd(4713.638236505974, 991.7657822506858));
  pnts.push_back(Vecd(4993.819530284301, 973.4675205855442));
  pnts.push_back(Vecd(5331.685208075813, 870.9972552607501));
  pnts.push_back(Vecd(5677.791512154923, 702.6532479414454));
  pnts.push_back(Vecd(6411.2072517511315, 574.5654162854526));
  pnts.push_back(Vecd(7408.323032550472, 464.77584629460216));
  pnts.push_back(Vecd(8298.310671611041, 358.64592863677944));
  pnts.push_back(Vecd(9163.576431808804, 300.0914913083252));
  pnts.push_back(Vecd(9666.254635352287, 142.72644098810588));
  pnts.push_back(Vecd(10012.360939431395, 300.0914913083261));
  pnts.push_back(Vecd(10234.857849196538, 387.9231473010068));
  pnts.push_back(Vecd(10564.482900700452, 402.5617566331198));
  pnts.push_back(Vecd(10688.09229501442, 483.0741079597433));
  pnts.push_back(Vecd(10861.145447053974, 753.8883806038425));
  pnts.push_back(Vecd(10877.62669962917, 563.5864592863677));
  pnts.push_back(Vecd(10762.2579316028, 336.6880146386093));
  pnts.push_back(Vecd(10358.467243510506, 190.3019213174748));
  pnts.push_back(Vecd(9855.789039967038, 62.21408966148147));
  pnts.push_back(Vecd(9715.698393077873, 25.61756633119785));
  pnts.push_back(Vecd(9460.23897816234, 29.277218664226893));
  pnts.push_back(Vecd(9130.613926658425, 117.10887465690757));
  pnts.push_back(Vecd(8751.545117428925, 190.3019213174748));
  pnts.push_back(Vecd(8281.829419035847, 204.94053064958825));
  pnts.push_back(Vecd(7573.135558302431, 263.4949679780416));
  pnts.push_back(Vecd(6716.110424392253, 365.96523330283617));
  pnts.push_back(Vecd(5718.994643592912, 468.4354986276303));
  pnts.push_back(Vecd(5166.872682323856, 629.4602012808782));
  pnts.push_back(Vecd(4730.119489081169, 768.526989935956));
  pnts.push_back(Vecd(3939.0193654717746, 834.4007319304665));
  pnts.push_back(Vecd(3156.1598681499786, 936.8709972552606));
  pnts.push_back(Vecd(2398.022249690976, 1134.4922232387921));
  pnts.push_back(Vecd(1747.0127729707456, 1427.2644098810615));
  pnts.push_back(Vecd(1153.6876802636998, 1632.2049405306493));
  pnts.push_back(Vecd(296.66254635352266, 2027.4473924977124));
  pnts.push_back(Vecd(8.240626287597706, 2371.4547118023784));
  pnts.push_back(Vecd(32.96250515039128, 3213.174748398902));
  return pnts;
}


class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWallShape(), ShapeBooleanOps::add);
    }
};
std::vector<Vecd> createSoilShape()
{
  std::vector<Vecd> pnts;
  pnts.push_back(Vecd(16.48125257519564, 3213.174748398902));
  pnts.push_back(Vecd(280.181293778327, 3103.3851784080507));
  pnts.push_back(Vecd(494.43757725587125, 2913.0832570905764));
  pnts.push_back(Vecd(766.3782447466006, 2686.184812442818));
  pnts.push_back(Vecd(1021.8376596621342, 2539.7987191216835));
  pnts.push_back(Vecd(1277.297074577668, 2305.5809698078683));
  pnts.push_back(Vecd(1532.7564894932011, 2107.959743824337));
  pnts.push_back(Vecd(1615.1627523691795, 1998.170173833486));
  pnts.push_back(Vecd(1747.0127729707456, 1851.784080512351));
  pnts.push_back(Vecd(1870.6221672847137, 1749.3138151875569));
  pnts.push_back(Vecd(1878.8627935723111, 1654.1628545288195));
  pnts.push_back(Vecd(1508.0346106304073, 1793.2296431838972));
  pnts.push_back(Vecd(1079.5220436753193, 1939.6157365050317));
  pnts.push_back(Vecd(576.8438401318499, 2225.068618481244));
  pnts.push_back(Vecd(271.9406674907291, 2612.9917657822507));
  pnts.push_back(Vecd(148.33127317676144, 2869.167429094236));
  pnts.push_back(Vecd(49.44375772558715, 3052.150045745654));
  pnts.push_back(Vecd(16.48125257519564, 3213.174748398902));
  return pnts;
}


class SoilBlock : public MultiPolygonShape
{
  public:
    explicit SoilBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createSoilShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createWallShape(), ShapeBooleanOps::sub);
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
#endif