function SolidStaticSchwarzController(params::Dict{Any,Any})
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["relative tolerance"]
    relative_tolerance = params["absolute tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    time = initial_time
    stop = 0
    prev_disp = Vector{Vector{Float64}}()
    SolidStaticSchwarzController(minimum_iterations, maximum_iterations, absolute_tolerance, relative_tolerance,
        initial_time, final_time, time_step, time, stop, prev_disp)
end

function SolidDynamicSchwarzController(params::Dict{Any,Any})
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["relative tolerance"]
    relative_tolerance = params["absolute tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    time = initial_time
    stop = 0
    prev_disp = Vector{Vector{Float64}}()
    prev_velo = Vector{Vector{Float64}}()
    prev_acce = Vector{Vector{Float64}}()
    SolidDynamicSchwarzController(minimum_iterations, maximum_iterations, absolute_tolerance, relative_tolerance,
        initial_time, final_time, time_step, time, stop, prev_disp, prev_velo, prev_acce)
end

function create_schwarz_controller(params::Dict{Any,Any})
    type = params["subdomains type"]
    if type == "static"
        return SolidStaticSchwarzController(params)
    elseif type == "dynamic"
        return SolidDynamicSchwarzController(params)
    else
        error("Unknown type of Schwarz controller : ", type)
    end
end

function schwarz(sim::MultiDomainSimulation)
    if sim.schwarz_controller.stop == 0
        solve(sim)
        write_step(sim)
    end
    set_subcycle_times(sim.schwarz_controller, sim.subsims)
    while true
        save_previous_solution(sim.schwarz_controller, sim.subsims)
        subcycle(sim.subsims)
    end
end

function save_previous_solution(schwarz_controller::SolidStaticSchwarzController, sims::Vector{SingleDomainSimulation})
    for sim ∈ sims
        push!(schwarz_controller.prev_disp, sim.integrator.displacement)
    end
end

function save_previous_solution(schwarz_controller::SolidDynamicSchwarzController, sims::Vector{SingleDomainSimulation})
    for sim ∈ sims
        push!(schwarz_controller.prev_disp, sim.integrator.displacement)
        push!(schwarz_controller.prev_velo, sim.integrator.velocity)
        push!(schwarz_controller.prev_acce, sim.integrator.acceleration)
    end
end

function set_subcycle_times(schwarz_controller::SchwarzController, sims::Vector{SingleDomainSimulation})
    initial_time = schwarz_controller.time
    final_time = round(schwarz_controller.time + schwarz_controller.time_step, digits=10)
    for sim ∈ sims
        sim.integrator.initial_time = initial_time
        sim.integrator.final_time = final_time
    end
end

function subcycle(sims::Vector{SingleDomainSimulation})
    for sim ∈ sims
        while continue_evolve(sim)
            apply_bcs(sim)
            advance(sim)
            advance_time(sim)
        end
    end
end
