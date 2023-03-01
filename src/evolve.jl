function evolve(simulation::Simulation)
    apply_ics(simulation)
    initialize_writing(simulation)
    while continue_evolve(simulation)
        synchronize(simulation)
        apply_bcs(simulation)
        advance(simulation)
        write_step(simulation)
        advance_time(simulation)
    end
    finalize_writing(simulation)
end

function continue_evolve(simulation::SingleDomainSimulation)
    return simulation.integrator.time ≤ simulation.integrator.final_time
end

function continue_evolve(simulation::MultiDomainSimulation)
    return simulation.schwarz_controller.time ≤ simulation.schwarz_controller.final_time
end

function advance(simulation::SingleDomainSimulation)
    solve(simulation)
end

function advance(simulation::MultiDomainSimulation)
    schwarz(simulation)
end

function apply_ics(simulation::SingleDomainSimulation)
    apply_ics(simulation.params, simulation.model)
end

function apply_ics(simulation::MultiDomainSimulation)
    for sub_simulation ∈ simulation.sub_simulations
        apply_ics(sub_simulation)
    end
end

function apply_bcs(simulation::SingleDomainSimulation)
    apply_bcs(simulation.params, simulation.model)
end

function apply_bcs(simulation::MultiDomainSimulation)
    for sub_simulation ∈ simulation.sub_simulations
        apply_bcs(sub_simulation)
    end
end

function solve(simulation::SingleDomainSimulation)
    solve(simulation.integrator, simulation.solver, simulation.model)
end

function solve(simulation::MultiDomainSimulation)
    for sub_simulation ∈ simulation.sub_simulations
        solve(sub_simulation)
    end
end

function initialize_writing(simulation::SingleDomainSimulation)
    initialize_writing(simulation.params, simulation.integrator, simulation.model)
end

function initialize_writing(simulation::MultiDomainSimulation)
    for sub_simulation ∈ simulation.sub_simulations
        initialize_writing(sub_simulation)
    end
end

function write_step(simulation::SingleDomainSimulation)
    println("Stop ", simulation.integrator.stop, ", time: ", simulation.integrator.time)
    write_step(simulation.params, simulation.integrator, simulation.model)
end

function write_step(simulation::MultiDomainSimulation)
    println("Stop ", simulation.schwarz_controller.stop, ", time: ", simulation.schwarz_controller.time)
    for sub_simulation ∈ simulation.sub_simulations
        write_step(sub_simulation)
    end
end

function finalize_writing(simulation::SingleDomainSimulation)
    finalize_writing(simulation.params)
end

function finalize_writing(simulation::MultiDomainSimulation)
    for sub_simulation ∈ simulation.sub_simulations
        finalize_writing(sub_simulation)
    end
end

function synchronize(simulation::SingleDomainSimulation)
    simulation.model.time = simulation.integrator.time
end

function synchronize(simulation::MultiDomainSimulation)
    time = simulation.schwarz_controller.time
    stop = simulation.schwarz_controller.stop
    for sub_simulation ∈ simulation.sub_simulations
        sub_simulation.integrator.time = sub_simulation.model.time = time
        sub_simulation.integrator.stop = stop
    end
end

function advance_time(simulation::SingleDomainSimulation)
    next_time = round(simulation.integrator.time + simulation.integrator.time_step; digits=10)
    simulation.integrator.time = next_time
    simulation.integrator.stop += 1
end

function advance_time(simulation::MultiDomainSimulation)
    next_time = round(simulation.schwarz_controller.time + simulation.schwarz_controller.time_step, digits=10)
    simulation.schwarz_controller.time = next_time
    simulation.schwarz_controller.stop += 1
end