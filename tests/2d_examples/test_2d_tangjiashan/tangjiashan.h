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
Real LL = 1500;                       /**< Soil column length. */
Real LH = 1500;                       /**< Soil column height. */
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
Real c_s = 2000; // sound speed
Real friction_angle = 30.0 * Pi / 180;
Real cohesion = 30e3;
Real dilatancy = 0.0;
 
//Diffusion
Real shear_period = 1.0;
std::string diffusion_species_name = "Concentration";
Real B = 1.0;
Real Sr = 0.016;
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
  pnts.push_back(Vecd(0.7038948850305076, 818.6536901865367));
  pnts.push_back(Vecd(34.490849366494615, 796.5936739659364));
  pnts.push_back(Vecd(107.69591740966686, 756.3665855636657));
  pnts.push_back(Vecd(178.08540591271705, 722.6277372262772));
  pnts.push_back(Vecd(252.69826372595028, 691.4841849148416));
  pnts.push_back(Vecd(317.4565931487565, 662.935928629359));
  pnts.push_back(Vecd(339.9812294697326, 652.5547445255472));
  pnts.push_back(Vecd(477.94462693571097, 649.9594484995943));
  pnts.push_back(Vecd(546.9263256687002, 684.9959448499592));
  pnts.push_back(Vecd(684.8897231346787, 784.9148418491484));
  pnts.push_back(Vecd(1086.109807602065, 1176.1557177615573));
  pnts.push_back(Vecd(1146.6447677146882, 1246.2287104622874));
  pnts.push_back(Vecd(1231.8160488033798, 1298.1346309813466));
  pnts.push_back(Vecd(1380.337869544815, 1336.4152473641527));
  pnts.push_back(Vecd(1442.2806194274992, 1337.7128953771294));
  pnts.push_back(Vecd(1432.426091037072, 1302.6763990267643));
  pnts.push_back(Vecd(1249.4134209291415, 1205.3527980535282));
  pnts.push_back(Vecd(1070.624120131394, 1085.969180859692));
  pnts.push_back(Vecd(908.7282965743783, 923.7631792376318));
  pnts.push_back(Vecd(766.5415297982169, 782.3195458231953));
  pnts.push_back(Vecd(631.3937118723604, 683.6982968369828));
  pnts.push_back(Vecd(508.9160018770531, 622.7088402270881));
  pnts.push_back(Vecd(448.3810417644298, 604.5417680454174));
  pnts.push_back(Vecd(348.4279680900986, 612.3276561232764));
  pnts.push_back(Vecd(218.9113092444862, 646.0665044606648));
  pnts.push_back(Vecd(114.73486625997188, 674.6147607461473));
  pnts.push_back(Vecd(28.85969028625061, 712.2465531224652));
  pnts.push_back(Vecd(0.7038948850305076, 727.8183292781829));
  pnts.push_back(Vecd(0.7038948850305076, 818.6536901865367));
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
  pnts.push_back(Vecd(414.5940872829658, 649.9594484995943));
  pnts.push_back(Vecd(448.3810417644298, 691.4841849148416));
  pnts.push_back(Vecd(493.43031440638197, 742.0924574209243));
  pnts.push_back(Vecd(572.2665415297984, 832.9278183292781));
  pnts.push_back(Vecd(624.3547630220554, 893.9172749391728));
  pnts.push_back(Vecd(649.6949788831536, 912.0843471208434));
  pnts.push_back(Vecd(682.0741435945566, 936.7396593673966));
  pnts.push_back(Vecd(731.3467855466919, 974.3714517437145));
  pnts.push_back(Vecd(810.1830126701079, 1031.4679643146796));
  pnts.push_back(Vecd(1017.1281088690756, 1175.506893755069));
  pnts.push_back(Vecd(1090.333176912248, 1213.138686131387));
  pnts.push_back(Vecd(1133.974659784139, 1233.9010543390107));
  pnts.push_back(Vecd(1136.7902393242612, 1214.4363341443636));
  pnts.push_back(Vecd(1101.5954950727362, 1157.3398215733982));
  pnts.push_back(Vecd(1055.138432660723, 1121.005677210057));
  pnts.push_back(Vecd(955.1853589863914, 1017.1938361719384));
  pnts.push_back(Vecd(848.1933364617552, 918.5725871857258));
  pnts.push_back(Vecd(720.0844673862038, 784.9148418491484));
  pnts.push_back(Vecd(637.0248709526045, 725.2230332522303));
  pnts.push_back(Vecd(548.3341154387613, 662.935928629359));
  pnts.push_back(Vecd(456.82778038479603, 639.5782643957823));
  pnts.push_back(Vecd(414.5940872829658, 649.9594484995943));
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