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
    apply_ics(simulation.sub_simulations)
    initialize_writing(simulation.sub_simulations)
    time = initial_time
    stop = 0
    time_step = params["time step"]
    while time <= final_time
        next_time = time + time_step
        synchronize(simulation.sub_simulations, stop, time)
        println("Stop ", stop, ", time: ", time)
        schwarz(simulation.sub_simulations, time, next_time, stop)
        time = round(next_time; digits=10)
        stop += 1
    end
    finalize_writing(models)
end

function schwarz(simulations::Vector{Simulation}, initial_time::Float64, final_time::Float64, stop::Int64)
    apply_bcs(simulations)
    if stop == 0
        solve(simulations)
        write_step(simulations)
    end
    write_step(simulations)
end

function apply_ics(simulation::Simulation)
    apply_ics(simulation.params, simulation.model)
end

function apply_bcs(simulation::Simulation)
    apply_bcs(simulation.params, simulation.model)
end

function apply_ics(simulations::Vector{Simulation})
    for simulation ∈ simulations
        apply_ics(simulation)
    end
end

function apply_bcs(simulations::Vector{Simulation})
    for simulation ∈ simulations
        apply_bcs(simulation)
    end
end

function solve(simulation::Simulation)
    solve(simulation.integrator, simulation.solver, simulation.model)
end

function solve(simulations::Vector{Simulation})
    for simulation ∈ simulations
        solve(simulation)
    end
end

function initialize_writing(simulation::Simulation)
    initialize_writing(simulation.params, simulation.integrator, simulation.model)
end

function initialize_writing(simulations::Vector{Simulation})
    for simulation ∈ simulations
        initialize_writing(simulation)
    end
end

function write_step(simulation::Simulation)
    write_step(simulation.params, simulation.integrator, simulation.model)
end

function write_step(simulations::Vector{Simulation})
    for simulation ∈ simulations
        write_step(simulation)
    end
end

function finalize_writing(simulation::Simulation)
    finalize_writing(simulation.params)
end

function finalize_writing(simulations::Vector{Simulation})
    for simulation ∈ simulations
        finalize_writing(simulation)
    end
end

function synchronize(simulations::Vector{Simulation}, time::Float64, stop::Int64)
    for simulation ∈ simulations
        simulation.integrator.time = simulation.model.time = time
        simulation.integrator.stop = stop
    end
end