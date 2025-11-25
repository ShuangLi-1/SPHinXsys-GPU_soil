/**
 * @file 	diffusion_NeumannBC.cpp
 * @brief 	2D test of diffusion problem with Neumann boundary condition.
 * @details This is the first case to validate multiple boundary conditions.
 * @author 	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */
#include "diffusion_DBC_CK.h"
#include "sphinxsys_ck.h"
using namespace SPH; // Namespace cite here
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<DiffusionBody>("DiffusionBody"));
    diffusion_body.defineMaterial<Solid>();
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_Dirichlet(sph_system, makeShared<DirichletWallBoundary>("DirichletWallBoundary"));
    wall_Dirichlet.defineMaterial<Solid>();
    wall_Dirichlet.generateParticles<BaseParticles, Lattice>();


    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    MultiPolygonShape lower_region_shape(MultiPolygon(lower_region_edge_points), "LowerRegion");
    BodyRegionByParticle wall_Dirichlet_lower_region(wall_Dirichlet, lower_region_shape);

    MultiPolygonShape upper_region_shape(MultiPolygon(upper_region_edge_points), "UpperRegion");
    BodyRegionByParticle wall_Dirichlet_upper_region(wall_Dirichlet, upper_region_shape);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Relation<Inner<>> diffusion_body_inner(diffusion_body);
    Relation<Contact<>> diffusion_body_contact_Dirichlet(diffusion_body, {&wall_Dirichlet});
    Relation<Contact<>> temperature_observer_contact(temperature_observer, {&diffusion_body});
    //----------------------------------------------------------------------
    // Define the main execution policy for this case.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy;
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> diffusion_body_cell_linked_list(diffusion_body);
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> wall_Dirichlet_cell_linked_list(wall_Dirichlet);

    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>>
        water_block_update_complex_relation(
            diffusion_body_inner, diffusion_body_contact_Dirichlet);
    UpdateRelation<MainExecutionPolicy, Contact<>> observer_contact_relation(temperature_observer_contact);

    StateDynamics<MainExecutionPolicy, InitialCondition<SPHBody, UniformDistribution<Real>>>
        diffusion_initial_condition(diffusion_body, diffusion_species_name, initial_temperature);
    StateDynamics<MainExecutionPolicy, InitialCondition<BodyRegionByParticle, UniformDistribution<Real>>>
        left_initial_condition(wall_Dirichlet_lower_region, diffusion_species_name, lower_temperature);
    StateDynamics<MainExecutionPolicy, InitialCondition<BodyRegionByParticle, UniformDistribution<Real>>>
        right_initial_condition(wall_Dirichlet_upper_region, diffusion_species_name, upper_temperature);

    IsotropicDiffusion isotropic_diffusion(diffusion_species_name, diffusion_coeff);
    GetDiffusionTimeStepSize get_time_step_size(diffusion_body, &isotropic_diffusion);
    RungeKuttaSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta1stStage, IsotropicDiffusion, KernelGradientInnerCK>,
            Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, KernelGradientContactCK>>,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta2ndStage, IsotropicDiffusion, KernelGradientInnerCK>,
            Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, KernelGradientContactCK>>>>
        diffusion_relaxation_rk2(DynamicsArgs(diffusion_body_inner, &isotropic_diffusion),
                                 DynamicsArgs(diffusion_body_contact_Dirichlet, &isotropic_diffusion));
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<MainExecutionPolicy, Real>>
        write_solid_temperature(diffusion_species_name, temperature_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------

    diffusion_body_cell_linked_list.exec();
    wall_Dirichlet_cell_linked_list.exec();

    water_block_update_complex_relation.exec();
    observer_contact_relation.exec();

    diffusion_initial_condition.exec();
    left_initial_condition.exec();
    right_initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 20;
    Real End_Time = T0;
    Real Observe_time = 0.005 * End_Time;
    Real Output_Time = 0.1 * End_Time;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TickCount::interval_t interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile(MainExecutionPolicy{});
    write_solid_temperature.writeToFile(ite);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < Output_Time)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observe_time)
            {
                if (ite % 500 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << sv_physical_time->getValue() << "	dt: "
                              << dt << "\n";
                }

                diffusion_relaxation_rk2.exec(dt);

                ite++;
                dt = get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                sv_physical_time->incrementValue(dt);
            }
            write_solid_temperature.writeToFile(ite);
        }

        TickCount t2 = TickCount::now();
        write_states.writeToFile(MainExecutionPolicy{});
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TickCount::interval_t tt;
    tt = t4 - t1 - interval;

    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Total physical time for computation: " << sv_physical_time->getValue() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_solid_temperature.generateDataBase(1.0e-3, 1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_solid_temperature.testResult();
    }

    return 0;
}