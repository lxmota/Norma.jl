function evolve(sim::Simulation)
    synchronize(sim)
    apply_ics(sim)
    apply_bcs(sim)
    initialize(sim)
    initialize_writing(sim)
    write_step(sim)
    while true
        advance_time(sim)
        if stop_evolve(sim) == true
            break
        end
        synchronize(sim)
        apply_bcs(sim)
        advance(sim)
        write_step(sim)
    end
    finalize_writing(sim)
end

function stop_evolve(sim::SingleDomainSimulation)
    return sim.integrator.time > sim.integrator.final_time
end

function stop_evolve(sim::MultiDomainSimulation)
    return sim.schwarz_controller.time > sim.schwarz_controller.final_time
end

function advance(sim::SingleDomainSimulation)
    solve(sim)
end

function advance(sim::MultiDomainSimulation)
    schwarz(sim)
end

function apply_ics(sim::SingleDomainSimulation)
    apply_ics(sim.params, sim.model)
end

function apply_ics(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        apply_ics(subsim)
    end
end

function apply_bcs(sim::SingleDomainSimulation)
    apply_bcs(sim.params, sim.model)
end

function apply_bcs(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        apply_bcs(subsim)
    end
end

function initialize(sim::SingleDomainSimulation)
    initialize(sim.integrator, sim.solver, sim.model)
end

function initialize(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        initialize(subsim)
    end
    type = sim.params["subdomains type"]
    if type == "dynamic solid mechanics"
        return
    elseif type == "static solid mechanics"
        schwarz(sim)
        return
    else
        error("Unknown type of Schwarz controller : ", type)
    end
end

function solve(sim::SingleDomainSimulation)
    solve(sim.integrator, sim.solver, sim.model)
end

function solve(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        solve(subsim)
    end
end

function initialize_writing(sim::SingleDomainSimulation)
    initialize_writing(sim.params, sim.integrator, sim.model)
end

function initialize_writing(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        initialize_writing(subsim)
    end
end

function write_step(sim::SingleDomainSimulation)
    write_step(sim.params, sim.integrator, sim.model)
end

function write_step(sim::MultiDomainSimulation)
    time = sim.schwarz_controller.time
    stop = sim.schwarz_controller.stop
    for subsim ∈ sim.subsims
        subsim.integrator.time = subsim.model.time = time
        subsim.integrator.stop = stop
        write_step(subsim)
    end
end

function finalize_writing(sim::SingleDomainSimulation)
    finalize_writing(sim.params)
end

function finalize_writing(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        finalize_writing(subsim)
    end
end

function synchronize(sim::SingleDomainSimulation)
    sim.model.time = sim.integrator.time
    println("stop ", sim.integrator.stop, ", t=", sim.integrator.time)
end

function synchronize(sim::MultiDomainSimulation)
    time = sim.schwarz_controller.time
    stop = sim.schwarz_controller.stop
    for subsim ∈ sim.subsims
        subsim.integrator.time = subsim.model.time = time
        subsim.integrator.stop = stop
    end
    println("stop ", sim.schwarz_controller.stop, ", t=", sim.schwarz_controller.time)
end

function advance_time(sim::SingleDomainSimulation)
    next_time = round(sim.integrator.time + sim.integrator.time_step; digits=10)
    sim.integrator.time = next_time
    sim.integrator.stop += 1
end

function advance_time(sim::MultiDomainSimulation)
    sim.schwarz_controller.prev_time = sim.schwarz_controller.time
    next_time = round(sim.schwarz_controller.time + sim.schwarz_controller.time_step, digits=10)
    sim.schwarz_controller.time = next_time
    sim.schwarz_controller.stop += 1
end