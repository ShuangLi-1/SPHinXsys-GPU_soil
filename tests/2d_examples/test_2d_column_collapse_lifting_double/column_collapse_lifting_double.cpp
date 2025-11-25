#include "column_collapse_lifting_double.h"

using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    std::cout << "diffusion_coeff:   " << diffusion_coeff <<std::endl;
    std::cout << "segregation_rate:   " << segregation_rate <<std::endl;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    RealBody soil_block(sph_system, makeShared<SoilBlock>("GranularBody"));
    //soil_block.defineClosure<PlasticContinuum, IsotropicDiffusion>(
        //ConstructArgs(rho0_s, c_s, Youngs_modulus, poisson, friction_angle,cohesion, dilatancy), ConstructArgs(diffusion_species_name, diffusion_coeff));
    soil_block.defineClosure<PlasticContinuum, SegregationLocalIsotropicDiffusion>(
        ConstructArgs(rho0_s, c_s, Youngs_modulus, poisson, friction_angle,cohesion, dilatancy, miu_s, miu_d), ConstructArgs(diffusion_species_name, diffusion_coeff, segregation_rate));
    soil_block.generateParticles<BaseParticles, Lattice>();
    
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation soil_block_inner(soil_block);
    ContactRelation soil_block_contact(soil_block, {&wall_boundary});
    ComplexRelation soil_block_complex(soil_block_inner, soil_block_contact);

    InnerRelation wall_boundary_inner(wall_boundary);
    ContactRelation wall_boundary_contact(wall_boundary, {&soil_block});
    ComplexRelation wall_boundary_complex(wall_boundary_inner, wall_boundary_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /*Geometric normal direction*/
    SimpleDynamics<NormalDirectionFromBodyShape> soil_block_normal_direction(soil_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /*Surface indicator*/
    InteractionWithUpdate<SegregationFreeSurfaceIndicationComplex> soil_block_indicator(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<SegregationFreeSurfaceIndicationComplex> wall_boundary_indicator(wall_boundary_inner, wall_boundary_contact);
     /*Plastic functions*/
    Real slope_angle = 0.0 * PI / 180;
    Gravity gravity(Vecd(gravity_g * sin(slope_angle), -gravity_g * cos(slope_angle)));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(soil_block, gravity);
    Dynamics1Level<continuum_dynamics::NonLocalPlasticIntegration1stHalfWithWallRiemann> granular_stress_relaxation(soil_block_inner, soil_block_contact);
    Dynamics1Level<continuum_dynamics::NonLocalPlasticIntegration2ndHalfWithWallRiemann> granular_density_relaxation(soil_block_inner, soil_block_contact);
    InteractionDynamics<continuum_dynamics::StressDiffusion> stress_diffusion(soil_block_inner);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> soil_acoustic_time_step(soil_block, 0.4);
    /*Update normal and degeneration*/
    Dynamics1Level<SoilSegAndPhiComplex> soil_update_seg_n_and_phi(soil_block_inner, soil_block_contact);
    Dynamics1Level<WallSegAndPhiComplex> wall_update_seg_n_and_phi(wall_boundary_inner, wall_boundary_contact);
    /*Kernel gradient correction*/
    //InteractionWithUpdate<LinearGradientCorrectionMatrixInner> soil_corrected_configuration(soil_block_inner);
    //InteractionWithUpdate<LinearGradientCorrectionMatrixInner> wall_corrected_configuration(wall_boundary_inner);
    
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> soil_corrected_configuration(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> wall_corrected_configuration(wall_boundary_inner, wall_boundary_contact);

    /*Diffusion*/
    CorrectedDiffusionBodyRelaxation temperature_relaxation(
        soil_block_inner, soil_block_contact);
    GetDiffusionTimeStepSize get_time_step_size(soil_block);
    SimpleDynamics<SheaboxDiffusionInitialCondition> setup_diffusion_initial_condition(soil_block);
    SimpleDynamics<NeumannWallBoundaryInitialCondition> setup_boundary_condition_Neumann(wall_boundary);
    
    /*Update diffusion coefficients*/
    SimpleDynamics<UpdateDiffusivityAndSegregationRate> update_D_and_V(soil_block, 1.0);
    SimpleDynamics<continuum_dynamics::UpdateNonLocalParameters> update_Non_Local_para(soil_block, 1.0);

    /*Lifting wall*/
    BodyRegionByParticle wall_lifting_maker(wall_boundary, makeShared<MultiPolygonShape>(createLiftingShape()));
    SimpleDynamics<Lifting> wall_lifting(wall_lifting_maker);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<int>(soil_block, "Indicator");
    write_states.addToWrite<Real>(soil_block, "SegPhi");
    write_states.addToWrite<Vecd>(soil_block, "SegNormalDirection");
    write_states.addToWrite<Matd>(soil_block, "LinearGradientCorrectionMatrix");
    write_states.addToWrite<Real>(soil_block, "ShearRate");
    write_states.addToWrite<Real>(soil_block, "Density");
    write_states.addToWrite<Real>(soil_block, "Pressure");
    write_states.addToWrite<Vecd>(soil_block, "Velocity");
    write_states.addToWrite<Real>(soil_block, "BPVisoosity");
    write_states.addToWrite<Real>(soil_block, "InertialNumber");
    write_states.addToWrite<Real>(soil_block, "AlphaPhi");
    write_states.addToWrite<Real>(soil_block, "Friction");
    write_states.addToWrite<Real>(soil_block, "Fluidity");
    write_states.addToWrite<Real>(soil_block, "LocalFluidityRate");
    write_states.addToWrite<Real>(soil_block, "NonLocalFluidityRate");

    write_states.addToWrite<int>(wall_boundary, "Indicator");
    write_states.addToWrite<Vecd>(wall_boundary, "SegNormalDirection");

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    setup_diffusion_initial_condition.exec();
    setup_boundary_condition_Neumann.exec();
    soil_block_normal_direction.exec();
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real End_Time =3.0;
    Real Output_Time = 0.05;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TickCount::interval_t interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < Output_Time)
        {
            Real relaxation_time = 0.0;

                if (ite % 500 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }
                stress_diffusion.exec();
                granular_stress_relaxation.exec(dt);
                granular_density_relaxation.exec(dt);

                soil_corrected_configuration.exec();
                wall_corrected_configuration.exec();

                soil_update_seg_n_and_phi.exec();

                soil_block_indicator.exec();
                wall_boundary_indicator.exec();

                update_Non_Local_para.exec();

                wall_lifting.exec();

                ite++;
                Real diffusion_dt = get_time_step_size.exec();
                Real dynamic_dt = soil_acoustic_time_step.exec();
                //dt = SMIN(diffusion_dt, dynamic_dt);
                dt = dynamic_dt;
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                soil_block.updateCellLinkedList();
                wall_boundary.updateCellLinkedList();
                soil_block_complex.updateConfiguration();
                wall_boundary_complex.updateConfiguration();
        }
        write_states.writeToFile();
        TickCount t2 = TickCount::now();

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TickCount::interval_t tt;
    tt = t4 - t1 - interval;

    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Total physical time for computation: " << physical_time << " seconds." << std::endl;

    return 0;
};
