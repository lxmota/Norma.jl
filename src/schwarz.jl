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

function schwarz(simulation::MultiDomainSimulation)
    if simulation.schwarz_controller.stop == 0
        solve(simulation)
        write_step(simulation)
    end
    set_subcycle_times(simulation.schwarz_controller, simulation.sub_simulations)
    while true
        save_previous_solution(simulation.schwarz_controller, simulation.sub_simulations)
        subcycle(simulation.sub_simulations)
    end
end

function save_previous_solution(schwarz_controller::SolidStaticSchwarzController, simulations::Vector{SingleDomainSimulation})
    for simulation ∈ simulations
        push!(schwarz_controller.prev_disp, simulation.integrator.displacement)
    end
end

function save_previous_solution(schwarz_controller::SolidDynamicSchwarzController, simulations::Vector{SingleDomainSimulation})
    for simulation ∈ simulations
        push!(schwarz_controller.prev_disp, simulation.integrator.displacement)
        push!(schwarz_controller.prev_velo, simulation.integrator.velocity)
        push!(schwarz_controller.prev_acce, simulation.integrator.acceleration)
    end
end

function set_subcycle_times(schwarz_controller::SchwarzController, simulations::Vector{SingleDomainSimulation})
    initial_time = schwarz_controller.time
    final_time = round(schwarz_controller.time + schwarz_controller.time_step, digits=10)
    for simulation ∈ simulations
        simulation.integrator.initial_time = initial_time
        simulation.integrator.final_time = final_time
    end
end

function subcycle(simulations::Vector{SingleDomainSimulation})
    for simulation ∈ simulations
        while continue_evolve(simulation)
            apply_bcs(simulation)
            advance(simulation)
            advance_time(simulation)
        end
    end
end
