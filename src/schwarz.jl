function SolidStaticSchwarzController(params::Dict{Any,Any})
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["relative tolerance"]
    relative_tolerance = params["absolute tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    absolute_error = 0.0
    relative_error = 0.0
    time = initial_time
    stop = 0
    converged = false
    prev_disp = Vector{Vector{Float64}}()
    SolidStaticSchwarzController(minimum_iterations, maximum_iterations, absolute_tolerance, relative_tolerance,
        absolute_error, relative_error, initial_time, final_time, time_step, time, stop, converged, prev_disp)
end

function SolidDynamicSchwarzController(params::Dict{Any,Any})
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["relative tolerance"]
    relative_tolerance = params["absolute tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    absolute_error = 0.0
    relative_error = 0.0
    time = initial_time
    stop = 0
    converged = false
    prev_disp = Vector{Vector{Float64}}()
    prev_velo = Vector{Vector{Float64}}()
    prev_acce = Vector{Vector{Float64}}()
    SolidDynamicSchwarzController(minimum_iterations, maximum_iterations, absolute_tolerance, relative_tolerance,
    absolute_error, relative_error, initial_time, final_time, time_step, time, stop, converged, prev_disp, prev_velo, prev_acce)
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
    set_subcycle_times(sim)
    iteration_number = 1
    while true
        save_previous_solutions(sim)
        subcycle(sim)
        iteration_number += 1
        update_convergence_criterion(sim)
        if continue_schwarz(sim, iteration_number) == false
            break
        end
    end
end

function save_previous_solutions(sim::MultiDomainSimulation)
    save_previous_solutions(sim.schwarz_controller, sim.subsims)
end

function save_previous_solutions(schwarz_controller::SolidStaticSchwarzController, sims::Vector{SingleDomainSimulation})
    for sim ∈ sims
        push!(schwarz_controller.prev_disp, sim.integrator.displacement)
    end
end

function save_previous_solutions(schwarz_controller::SolidDynamicSchwarzController, sims::Vector{SingleDomainSimulation})
    for sim ∈ sims
        push!(schwarz_controller.prev_disp, sim.integrator.displacement)
        push!(schwarz_controller.prev_velo, sim.integrator.velocity)
        push!(schwarz_controller.prev_acce, sim.integrator.acceleration)
    end
end

function set_subcycle_times(sim::MultiDomainSimulation)
    initial_time = sim.schwarz_controller.time
    final_time = round(sim.schwarz_controller.time + sim.schwarz_controller.time_step, digits=10)
    for subsim ∈ sim.subsims
        subsim.integrator.initial_time = initial_time
        subsim.integrator.final_time = final_time
    end
end

function subcycle(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        while continue_evolve(subsim)
            apply_bcs(subsim)
            advance(subsim)
            advance_time(subsim)
        end
    end
end

function update_convergence_criterion(sim::MultiDomainSimulation)
    sim.schwarz_controller.converged = is_converged(sim.schwarz_controller, sim.subsims)
end

function is_converged(schwarz_controller::SolidStaticSchwarzController, sims::Vector{SingleDomainSimulation})
    num_sims = length(sims)
    norms_disp = zeros(num_sims)
    norms_diff = zeros(num_sims)
    for i ∈ 1:num_sims
        xᵖʳᵉᵛ = schwarz_controller.prev_disp[i]
        xᵖʳᵉᵛ = sims[i].integrator.displacement
        norms_disp[i] = norm(xᵖʳᵉᵛ)
        norms_diff[i] = norm(xᶜᵘʳʳ - xᵖʳᵉᵛ)
    end
    norm_disp = norm(norms_disp)
    norm_diff = norm(norms_diff)
    schwarz_controller.absolute_error = norm_diff
    schwarz_controller.relative_error = norm_disp > 0.0 ? norm_diff / norm_disp : norm_diff
    conv_abs = schwarz_controller.absolute_error ≤ schwarz_controller.absolute_tolerance
    conv_rel = schwarz_controller.relative_error ≤ schwarz_controller.relative_tolerance
    return conv_abs || conv_rel
end

function is_converged(schwarz_controller::SolidDynamicSchwarzController, sims::Vector{SingleDomainSimulation})
    num_sims = length(sims)
    norms_disp = zeros(num_sims)
    norms_diff = zeros(num_sims)
    for i ∈ 1:num_sims
        Δt = schwarz_controller.time_step
        xᵖʳᵉᵛ = schwarz_controller.prev_disp[i] + Δt * schwarz_controller.prev_velo[i]
        xᵖʳᵉᵛ = sims[i].integrator.displacement + Δt * sims[i].integrator.velocity
        norms_disp[i] = norm(xᵖʳᵉᵛ)
        norms_diff[i] = norm(xᵖʳᵉᵛ - xᵖʳᵉᵛ)
    end
    norm_disp = norm(norms_disp)
    norm_diff = norm(norms_diff)
    schwarz_controller.absolute_error = norm_diff
    schwarz_controller.relative_error = norm_disp > 0.0 ? norm_diff / norm_disp : norm_diff
    conv_abs = schwarz_controller.absolute_error ≤ schwarz_controller.absolute_tolerance
    conv_rel = schwarz_controller.relative_error ≤ schwarz_controller.relative_tolerance
    return conv_abs || conv_rel
end

function continue_schwarz(sim::MultiDomainSimulation, iteration_number::Int64)
    for subsim ∈ sim.subsims
        if subsim.solver.failed == true
            return false
        end
    end
    if sim.schwarz_controller.absolute_error == 0.0
        return false
    end
    exceeds_minimum_iterations = iteration_number > sim.schwarz_controller.minimum_iterations
    if exceeds_minimum_iterations == false
        return true
    end
    exceeds_maximum_iterations = iteration_number > sim.schwarz_controller.maximum_iterations
    if exceeds_maximum_iterations == true
        return false
    end
    unconverged = sim.schwarz_controller.converged == false
    return unconverged
end