/**
 * @file 	column_collapse.cpp
 * @brief 	3D repose angle example.
 * @details This is the one of the basic test cases, also the first case for understanding
 * 			SPH method for modelling granular materials such as soils and sands.
 * @author Shuaihao Zhang and Xiangyu Hu
 */
#include "damFormation_zhou_non_local.h"
using namespace SPH;


int main(int ac, char *av[])
{
    //----------------------------------------------------------------------shear
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vecd(-0.5 * DL - BW, -DH - BW, -DW - BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref, 128);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();

    std::cout<<"cs:  " <<c_s<<std::endl;
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    RealBody soil_block(sph_system, makeShared<SoilBlock>("GranularBody"));
    soil_block.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    soil_block.defineClosure<PlasticContinuum, SegregationLocalIsotropicDiffusion>(
        ConstructArgs(rho0_s, c_s, Youngs_modulus, poisson, friction_angle,cohesion, dilatancy, miu_s, miu_d, d_min ,d_max), ConstructArgs(diffusion_species_name, diffusion_coeff, segregation_rate));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? soil_block.generateParticles<BaseParticles, Reload>(soil_block.getName())
        : soil_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    wall_boundary.defineMaterial<Solid>();
    if(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {     
        wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName());
        std::cout<<"reload wall boundary"<<std::endl;
    }
    else
    {     wall_boundary.generateParticles<BaseParticles, Lattice>();}
     //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation wall_boundary_inner(wall_boundary);
        InnerRelation soil_block_inner(soil_block);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> wall_random_particles(wall_boundary);
        SimpleDynamics<RandomizeParticlePosition> soil_random_particles(soil_block);
        BodyStatesRecordingToVtp write_body_to_vtp(sph_system);
        ReloadParticleIO write_particle_reload_files({&wall_boundary, &soil_block});
        RelaxationStepInner wall_relaxation_step_inner(wall_boundary_inner);
        RelaxationStepInner soil_relaxation_step_inner(soil_block_inner);

        // ----------------------------------------------------------------------
        // 	Particle relaxation starts here.
        // ----------------------------------------------------------------------
        wall_random_particles.exec(0.25);
        soil_random_particles.exec(0.25);
        wall_relaxation_step_inner.SurfaceBounding().exec();
        soil_relaxation_step_inner.SurfaceBounding().exec();
        write_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            wall_relaxation_step_inner.exec();
            soil_relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation soil_block_inner(soil_block);
    ContactRelation soil_block_contact(soil_block, {&wall_boundary});
    ComplexRelation soil_block_complex(soil_block_inner, soil_block_contact);

    InnerRelation wall_boundary_inner(wall_boundary);
    ContactRelation wall_boundary_contact(wall_boundary, {&soil_block});
    ComplexRelation wall_boundary_complex(wall_boundary_inner, wall_boundary_contact);
     /*Geometric normal direction*/
    SimpleDynamics<NormalDirectionFromBodyShape> soil_block_normal_direction(soil_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /*Surface indicator*/
    InteractionWithUpdate<SegregationFreeSurfaceIndicationComplex> soil_block_indicator(soil_block_inner, soil_block_contact);
    //InteractionWithUpdate<SegregationFreeSurfaceIndicationComplex> wall_boundary_indicator(wall_boundary_inner, wall_boundary_contact);
    InteractionWithUpdate<FreeSurfaceIndicationComplex> wall_boundary_indicator(wall_boundary_inner, wall_boundary_contact);
     /*Plastic functions*/
    Real slope_angle = 0.0 * PI / 180;
    Gravity gravity(Vecd(gravity_g * sin(slope_angle), -gravity_g * cos(slope_angle),0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(soil_block, gravity);
    Dynamics1Level<continuum_dynamics::NonLocalPlasticIntegration1stHalfWithWallRiemann> granular_stress_relaxation(soil_block_inner, soil_block_contact);
    Dynamics1Level<continuum_dynamics::NonLocalPlasticIntegration2ndHalfWithWallRiemann> granular_density_relaxation(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> soil_density_by_summation(soil_block_inner, soil_block_contact);
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
    DiffusionBodyRelaxation temperature_relaxation(
        soil_block_inner, soil_block_contact);
    GetDiffusionTimeStepSize get_time_step_size(soil_block);
    SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(soil_block);
    SimpleDynamics<NeumannWallBoundaryInitialCondition> setup_boundary_condition_Neumann(wall_boundary);
    
    /*Update diffusion coefficients*/
    SimpleDynamics<UpdateDiffusivityAndSegregationRate> update_D_and_V(soil_block, 1.0);
    SimpleDynamics<continuum_dynamics::UpdateNonLocalParameters> update_Non_Local_para(soil_block, 1.0);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<int>(soil_block, "Indicator");
    write_states.addToWrite<Real>(soil_block, "Test");
    write_states.addToWrite<Real>(soil_block, "SegPhi");
    write_states.addToWrite<Vecd>(soil_block, "SegNormalDirection");
    write_states.addToWrite<Matd>(soil_block, "LinearGradientCorrectionMatrix");
    write_states.addToWrite<Real>(soil_block, "ShearRate");
    write_states.addToWrite<Real>(soil_block, "Density");
    write_states.addToWrite<Vecd>(soil_block, "Velocity");
    write_states.addToWrite<Real>(soil_block, "Pressure");
    write_states.addToWrite<Real>(soil_block, "InertialNumber");
    write_states.addToWrite<Real>(soil_block, "Friction");
    write_states.addToWrite<Real>(soil_block, "Fluidity");
    write_states.addToWrite<Real>(soil_block, "LocalFluidityRate");
    write_states.addToWrite<Real>(soil_block, "NonLocalFluidityRate");
    write_states.addToWrite<Real>(soil_block, "AlphaPhi");
    write_states.addToWrite<Real>(soil_block, "KC");
    write_states.addToWrite<Real>(soil_block, "ConcentrationChangeRate");

    write_states.addToWrite<int>(wall_boundary, "Indicator");
    write_states.addToWrite<Vecd>(wall_boundary, "SegNormalDirection");
    write_states.addToWrite<Vecd>(wall_boundary, "NormalDirection");

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

    /*Wall update only once*/
    wall_corrected_configuration.exec();
    wall_boundary_indicator.exec();
    wall_update_seg_n_and_phi.exec();
    update_D_and_V.exec();
    update_Non_Local_para.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 150 * shear_period;
    Real End_Time = 3.0;
    Real Output_Time = 0.1;
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
                soil_density_by_summation.exec();
                stress_diffusion.exec();
                granular_stress_relaxation.exec(dt);
                granular_density_relaxation.exec(dt);

                // soil_corrected_configuration.exec();
                // //wall_corrected_configuration.exec();


                // soil_block_indicator.exec();
                // //wall_boundary_indicator.exec();

                // soil_update_seg_n_and_phi.exec();
                // //wall_update_seg_n_and_phi.exec();

                // update_D_and_V.exec();
                // update_Non_Local_para.exec();

                // temperature_relaxation.exec(dt);

                ite++;
                Real diffusion_dt = get_time_step_size.exec();
                Real dynamic_dt = soil_acoustic_time_step.exec();
                dt = SMIN(diffusion_dt, dynamic_dt);
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                soil_block.updateCellLinkedList();
                //wall_boundary.updateCellLinkedList();
                soil_block_complex.updateConfiguration();
                //wall_boundary_complex.updateConfiguration();
        }
        // write_states.writeToFile();
        TickCount t2 = TickCount::now();

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TickCount::interval_t tt;
    tt = t4 - t1 - interval;
    write_states.writeToFile();
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Total physical time for computation: " << physical_time << " seconds." << std::endl;

    return 0;
}
