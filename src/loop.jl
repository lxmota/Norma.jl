function loop(simulation::SingleDomainSimulation)
    apply_ics(simulation)
    initialize_writing(simulation)
    while simulation.integrator.time <= simulation.integrator.final_time
        simulation.model.time = simulation.integrator.time
        println("Stop ", simulation.integrator.stop, ", time: ", simulation.integrator.time)
        apply_bcs(simulation)
        solve(simulation)
        write_step(simulation)
        simulation.integrator.time = round(simulation.integrator.time + simulation.integrator.time_step; digits=10)
        simulation.integrator.stop += 1
    end
    finalize_writing(simulation)
end

function loop(simulation::MultiDomainSimulation)
    apply_ics(simulation)
    initialize_writing(simulation)
    final_time = simulation.schwarz_controller.final_time
    time_step = simulation.schwarz_controller.time_step
    while simulation.schwarz_controller.time <= final_time
        next_time = time + time_step
        synchronize(simulation)
        println("Stop ", stop, ", time: ", time)
        schwarz(simulation, time, next_time, stop)
        simulation.schwarz_controller.time = round(next_time; digits=10)
        stop += 1
    end
    finalize_writing(simulation)
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
    write_step(simulation.params, simulation.integrator, simulation.model)
end

function write_step(simulation::MultiDomainSimulation)
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

function synchronize(simulation::MultiDomainSimulation)
    time = simulation.schwarz_controller.time
    stop = simulation.schwarz_controller.stop
    for sub_simulation ∈ simulation.sub_simulations
        sub_simulation.integrator.time = sub_simulation.model.time = time
        sub_simulation.integrator.stop = stop
    end
end