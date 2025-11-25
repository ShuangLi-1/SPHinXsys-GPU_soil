/**
 * @file 	segregation_ck.cpp
 * @brief 	2D particle segregation.
 * @details The method from Zhu et al 2023 (JFM)
 * @author Shuang Li, Xiangyu Hu 
 */
#include "sphinxsys_ck.h" //SPHinXsys Library.
using namespace SPH;
#include "test_model.h" //header for this case
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
    soil_block.defineClosure<PlasticContinuum, LocalIsotropicDiffusion>(
        ConstructArgs(rho0_s, c_s, Youngs_modulus, poisson, friction_angle,cohesion, dilatancy), ConstructArgs(diffusion_species_name, diffusion_coeff, segregation_rate));
    soil_block.generateParticles<BaseParticles, Lattice>();
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    
    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Creating body parts.
    //---------------------------------------------------------------------
    MultiPolygonShape lower_region_shape(MultiPolygon(lower_region_edge_points), "LowerRegion");
    BodyRegionByParticle wall_Dirichlet_lower_region(wall_boundary, lower_region_shape);
    MultiPolygonShape upper_region_shape(MultiPolygon(upper_region_edge_points), "UpperRegion");
    BodyRegionByParticle wall_Dirichlet_upper_region(wall_boundary, upper_region_shape);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy; // define execution policy for this case

    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> soil_cell_linked_list(soil_block);
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> wall_cell_linked_list(wall_boundary);

    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    //No ck interaction
    InnerRelation wall_boundary_inner_old(wall_boundary);
    InnerRelation soil_inner_old(soil_block);

    Relation<Inner<>> soil_block_inner(soil_block);
    Relation<Contact<>> soil_block_contact(soil_block, {&wall_boundary});
    Relation<Inner<>>  wall_boundary_inner(wall_boundary);
    Relation<Contact<>> wall_boundary_contact(wall_boundary, {&soil_block});
    Relation<Contact<>> temperature_observer_contact(temperature_observer, {&soil_block});

    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> soil_block_update_complex_relation(soil_block_inner, soil_block_contact);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> wall_boundary_update_complex_relation(wall_boundary_inner, wall_boundary_contact);
        UpdateRelation<MainExecutionPolicy, Contact<>> observer_contact_relation(temperature_observer_contact);
    ParticleSortCK<MainExecutionPolicy, QuickSort> particle_sort(soil_block);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    StateDynamics<MainExecutionPolicy, GravityForceAndPos0CK<Gravity>> constant_gravity(soil_block, gravity);
    //StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<NormalDirectionFromParticles> wall_boundary_normal_direction(wall_boundary_inner_old);
    InteractionDynamics<NormalDirectionFromParticles> soil_normal_direction(soil_inner_old);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> soil_advection_step_setup(soil_block);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepClose> soil_advection_step_close(soil_block);

    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::PlasticAcousticStep1stHalfWithWallRiemannCK>
         soil_acoustic_step_1st_half(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::PlasticAcousticStep2ndHalfWithWallRiemannCK>
        soil_acoustic_step_2nd_half(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensityRegularizationComplexFreeSurface>
        soil_density_regularization(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::StressDiffusionInnerCK> stress_diffusion(soil_block_inner);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::SegFreeSurfaceIndicationComplexCK>
        soil_boundary_indicator(soil_block_inner, soil_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::SegFreeSurfaceIndicationComplexCK>
        wall_boundary_indicator(wall_boundary_inner, wall_boundary_contact);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK> soil_acoustic_time_step(soil_block,0.4);
    //StateDynamics<MainExecutionPolicy, fluid_dynamics::ShearboxMovement> soil_shear_box_movement(soil_block);
    /*Kernel Correction*/
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
    soil_linear_correction_matrix(DynamicsArgs(soil_block_inner, 0.5), soil_block_contact);
    //----------------------------------------------------------------------//
    // Diffusion functions
    //----------------------------------------------------------------------//
    //Set initilization
    StateDynamics<MainExecutionPolicy, InitialCondition<SPHBody, UniformDistribution<Real>>>
        diffusion_initial_condition(soil_block, diffusion_species_name, initial_temperature);
    StateDynamics<MainExecutionPolicy, InitialCondition<BodyRegionByParticle, UniformDistribution<Real>>>
        lower_initial_condition(wall_Dirichlet_lower_region, diffusion_species_name, lower_temperature);
    StateDynamics<MainExecutionPolicy, InitialCondition<BodyRegionByParticle, UniformDistribution<Real>>>
        upper_initial_condition(wall_Dirichlet_upper_region, diffusion_species_name, upper_temperature);
    //StateDynamics<MainExecutionPolicy, continuum_dynamics::SegregationParametersCK> SegregationParameter(soil_block);
    
    //Dynamics
    LocalIsotropicDiffusion isotropic_diffusion(diffusion_species_name, diffusion_coeff, diffusion_coeff);
    GetDiffusionTimeStepSize get_time_step_size(soil_block, &isotropic_diffusion);
    /*Segregation Normal Direction*/
    // InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::UpdateSegNormalDirectionRiemann>
    //     update_seg_n(soil_block_inner, soil_block_contact);
    // InteractionDynamicsCK<MainExecutionPolicy, continuum_dynamics::WallUpdateSegNormalDirectionRiemann>
    //     wall_update_seg_n(wall_boundary_inner, wall_boundary_contact);
    /*Main diffusion function*/
    RungeKuttaSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta1stStage, LocalIsotropicDiffusion, KernelGradientInnerCK>,
            Contact<InteractionOnly, Dirichlet<LocalIsotropicDiffusion>, KernelGradientContactCK>>,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta2ndStage, LocalIsotropicDiffusion, KernelGradientInnerCK>,
            Contact<InteractionOnly, Dirichlet<LocalIsotropicDiffusion>, KernelGradientContactCK>>>>
        diffusion_relaxation_rk2(DynamicsArgs(soil_block_inner, &isotropic_diffusion),
                                 DynamicsArgs(soil_block_contact, &isotropic_diffusion));
    //----------------------------------------------------------------------
    //	Define shear motion.
    //----------------------------------------------------------------------
    //BodyRegionByParticle shear_maker(wall_boundary, makeShared<MultiPolygonShape>(createShearShape()));
    //SimpleDynamics<ShearMaking> shear_making(shear_maker);
    //StateDynamics<MainExecutionPolicy, ShearMaking> shear_making(shear_maker);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    //body_states_recording.addToWrite<Vecd>(wall_boundary, "SegNormalDirection");

    body_states_recording.addToWrite<Real>(soil_block, "Density");
    StateDynamics<MainExecutionPolicy,continuum_dynamics::VerticalStressCK> vertical_stress(soil_block);
    body_states_recording.addToWrite<Real>(soil_block, "VerticalStress");
    StateDynamics<MainExecutionPolicy,continuum_dynamics::AccDeviatoricPlasticStrainCK> accumulated_deviatoric_plastic_strain(soil_block);
    body_states_recording.addToWrite<Real>(soil_block, "AccDeviatoricPlasticStrain");
    body_states_recording.addToWrite<Vecd>(soil_block, "InitialPosition");

    //body_states_recording.addToWrite<Vecd>(soil_block, "NormalDirection");
    // body_states_recording.addToWrite<Real>(soil_block, "SegregationDiffusivity");
    // body_states_recording.addToWrite<Vecd>(soil_block, "SegNormalDirection");
    // body_states_recording.addToWrite<int>(soil_block, "Indicator");
    
    RestartIO restart_io(sph_system);
    
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<MainExecutionPolicy, TotalMechanicalEnergyCK>>
    write_mechanical_energy(soil_block, gravity);

    RegressionTestEnsembleAverage<ObservedQuantityRecording<MainExecutionPolicy, Real>>
    write_soil_temperature(diffusion_species_name, temperature_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");    

    constant_gravity.exec();
    soil_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    soil_block_update_complex_relation.exec();
    wall_boundary_update_complex_relation.exec();
    observer_contact_relation.exec();

    sph_system.initializeSystemCellLinkedLists();
    // initial periodic boundary condition
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    soil_normal_direction.exec();
    //Diffusion initilization
    diffusion_initial_condition.exec();
    lower_initial_condition.exec();
    upper_initial_condition.exec();
    //wall_Neumann_initial_condition.exec();
    //SegregationParameter.exec();

    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------    
    size_t number_of_iterations = 0;
    int screen_output_interval = 10;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real End_Time = 20.0;         /**< End time. */
    Real D_Time = 0.01; /**< Time stamps for output of body states. */
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
    write_soil_temperature.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue()  < End_Time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            //shear_making.exec();
            /** outer loop for dual-time criteria time-stepping. */
            soil_density_regularization.exec();
            soil_boundary_indicator.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
                    
            //soil_advection_step_setup.exec();
            Real dt = get_time_step_size.exec();
            //stress_diffusion.exec();
            //soil_acoustic_step_1st_half.exec(dt);
            //soil_acoustic_step_2nd_half.exec(dt);
            //soil_shear_box_movement.exec(dt);
            // update_seg_n.exec(dt);
            // wall_update_seg_n.exec(dt);
            diffusion_relaxation_rk2.exec(dt);
            //soil_advection_step_close.exec();
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
                    write_soil_temperature.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            
            number_of_iterations++;
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();              
            soil_cell_linked_list.exec();
            wall_cell_linked_list.exec();
            soil_block_update_complex_relation.exec();
            wall_boundary_update_complex_relation.exec();
            observer_contact_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            wall_boundary.updateCellLinkedList();
            wall_boundary_inner_old.updateConfiguration();
            soil_inner_old.updateConfiguration();
            wall_boundary_normal_direction.exec();
            soil_normal_direction.exec();
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
