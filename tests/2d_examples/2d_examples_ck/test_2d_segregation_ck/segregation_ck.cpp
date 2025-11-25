/**
 * @file 	segregation_ck.cpp
 * @brief 	2D particle segregation.
 * @details The method from Zhu et al 2023 (JFM)
 * @author Shuang Li, Xiangyu Hu 
 */
#include "sphinxsys_ck.h" //SPHinXsys Library.
using namespace SPH;
#include "segregation_ck.h" //header for this case
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
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
        ConstructArgs(rho0_s, c_s, Youngs_modulus, poisson, friction_angle,cohesion, dilatancy), ConstructArgs(diffusion_species_name, diffusion_coeff, segregation_rate));
    soil_block.generateParticles<BaseParticles, Lattice>();

    FluidBody water_block(sph_system, makeShared<SoilBlock>("FluidBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, SegregationLocalIsotropicDiffusion>(
        ConstructArgs(rho0_s,c_s), ConstructArgs(diffusion_species_name, diffusion_coeff, segregation_rate));
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy; // define execution policy for this case

    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> soil_cell_linked_list(soil_block);
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> wall_cell_linked_list(wall_boundary);

    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    Relation<Inner<>> soil_block_inner(soil_block);
    Relation<Contact<>> soil_block_contact(soil_block, {&wall_boundary});

    Relation<Inner<>> water_block_inner(water_block);
    Relation<Contact<>> water_block_contact(water_block, {&wall_boundary});

    InnerRelation wall_boundary_inner(wall_boundary);

    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> soil_block_update_complex_relation(soil_block_inner, soil_block_contact);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_block_update_complex_relation(water_block_inner, water_block_contact);
    ParticleSortCK<MainExecutionPolicy, QuickSort> particle_sort(soil_block);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    StateDynamics<MainExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(soil_block, gravity);
    //StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<NormalDirectionFromParticles> wall_boundary_normal_direction(wall_boundary_inner);
    /*Soil dynamics*/
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> soil_advection_step_setup(soil_block);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepClose> soil_advection_step_close(soil_block);
    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::PlasticAcousticStep1stHalfWithWallRiemannCK>
        soil_acoustic_step_1st_half(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::PlasticAcousticStep2ndHalfWithWallRiemannCK>
        soil_acoustic_step_2nd_half(soil_block_inner, soil_block_contact);
    /*Fluid dynamics*/
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_block);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepClose> water_advection_step_close(water_block);
        InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(water_block_inner, 0.5), water_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_1st_half(water_block_inner, water_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_2nd_half(water_block_inner, water_block_contact);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensityRegularizationComplexFreeSurface>
        soil_density_regularization(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::StressDiffusionInnerCK> stress_diffusion(soil_block_inner);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK> soil_acoustic_time_step(soil_block,0.4);
    
    //----------------------------------------------------------------------//
    // Diffusion functions
    //----------------------------------------------------------------------//
    //Set initilization
    StateDynamics<MainExecutionPolicy, InitialCondition<SPHBody, ShearBoxDistribution<Real>>>
        diffusion_initial_condition(soil_block, diffusion_species_name, initial_temperature);
    StateDynamics<MainExecutionPolicy, InitialCondition<SPHBody, UniformDistribution<Real>>>
        wall_Neumann_initial_condition(wall_boundary, diffusion_species_name + "Flux", heat_flux);
    StateDynamics<MainExecutionPolicy, continuum_dynamics::SegregationParametersCK> SegregationParameter(soil_block);
    
    //Dynamics
    SegregationLocalIsotropicDiffusion isotropic_diffusion(diffusion_species_name, diffusion_coeff, diffusion_coeff);
    GetDiffusionTimeStepSize get_time_step_size(soil_block, &isotropic_diffusion);
    RungeKuttaSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        SegregationDiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta1stStage, SegregationLocalIsotropicDiffusion, KernelGradientInnerCK>,
            Contact<InteractionOnly, Neumann<SegregationLocalIsotropicDiffusion>, KernelGradientContactCK>>,
        SegregationDiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta2ndStage, SegregationLocalIsotropicDiffusion, KernelGradientInnerCK>,
            Contact<InteractionOnly, Neumann<SegregationLocalIsotropicDiffusion>, KernelGradientContactCK>>>>
        diffusion_relaxation_rk2(DynamicsArgs(soil_block_inner, &isotropic_diffusion),
                                 DynamicsArgs(soil_block_contact, &isotropic_diffusion));
    //----------------------------------------------------------------------
    //	Define shear motion.
    //----------------------------------------------------------------------
    BodyRegionByParticle shear_maker(wall_boundary, makeShared<MultiPolygonShape>(createShearShape()));
    SimpleDynamics<ShearMaking> shear_making(shear_maker);
    //StateDynamics<MainExecutionPolicy, ShearMaking> shear_making(shear_maker);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Real>(soil_block, "Density");
    StateDynamics<MainExecutionPolicy,continuum_dynamics::VerticalStressCK> vertical_stress(soil_block);
    body_states_recording.addToWrite<Real>(soil_block, "VerticalStress");
    StateDynamics<MainExecutionPolicy,continuum_dynamics::AccDeviatoricPlasticStrainCK> accumulated_deviatoric_plastic_strain(soil_block);
    body_states_recording.addToWrite<Real>(soil_block, "AccDeviatoricPlasticStrain");

    body_states_recording.addToWrite<Real>(soil_block, "SegregationDiffusivity");
    
    RestartIO restart_io(sph_system);
    
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<MainExecutionPolicy, TotalMechanicalEnergyCK>>
    write_mechanical_energy(soil_block, gravity);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");    

    constant_gravity.exec();
    soil_cell_linked_list.exec();
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    soil_block_update_complex_relation.exec();
    water_block_update_complex_relation.exec();

    sph_system.initializeSystemCellLinkedLists();
    // initial periodic boundary condition
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    //Diffusion initilization
    diffusion_initial_condition.exec();
    wall_Neumann_initial_condition.exec();
    //SegregationParameter.exec();

    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------    
    size_t number_of_iterations = 0;
    int screen_output_interval = 500;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real End_Time = 1300.0;         /**< End time. */
    Real D_Time = 1.0; /**< Time stamps for output of body states. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_acoustic_steps;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(MainExecutionPolicy{});
    write_mechanical_energy.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue()  < End_Time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            shear_making.exec();
            /** outer loop for dual-time criteria time-stepping. */
            soil_density_regularization.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
                    

            Real dt = soil_acoustic_time_step.exec();

            /*Soil dynamics*/            
            // soil_advection_step_setup.exec();
            // stress_diffusion.exec();
            // soil_acoustic_step_1st_half.exec(dt);
            // soil_acoustic_step_2nd_half.exec(dt);
            // soil_advection_step_close.exec();
            /*Fluid dynamics*/
            water_advection_step_setup.exec();
            fluid_linear_correction_matrix.exec();
            fluid_acoustic_step_1st_half.exec(dt);
            fluid_acoustic_step_2nd_half.exec(dt);
            water_advection_step_close.exec();

            /*Diffusion model*/
            diffusion_relaxation_rk2.exec(dt);
            relaxation_time += dt;
            integration_time += dt;
            sv_physical_time->incrementValue(dt);

            interval_acoustic_steps += TickCount::now() - time_instance;

            /** screen output, write body observables and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << std::setprecision(4) << "	Time = "
                            << sv_physical_time->getValue()
                            << std::scientific << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_mechanical_energy.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();              
            soil_cell_linked_list.exec();
            water_cell_linked_list.exec();
            soil_block_update_complex_relation.exec();
            water_block_update_complex_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            wall_boundary.updateCellLinkedList();
            wall_boundary_inner.updateConfiguration();
            wall_boundary_normal_direction.exec();
        }
        vertical_stress.exec();
        accumulated_deviatoric_plastic_strain.exec();
        body_states_recording.writeToFile(MainExecutionPolicy{});
        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
            
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << std::fixed << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_mechanical_energy.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_mechanical_energy.testResult();
    }

    return 0;
};
