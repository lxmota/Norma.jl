function SolidSchwarzController(params::Dict{Any,Any})
    num_domains = length(params["domains"])
    minimum_iterations = params["minimum iterations"]
    maximum_iterations = params["maximum iterations"]
    absolute_tolerance = params["relative tolerance"]
    relative_tolerance = params["absolute tolerance"]
    initial_time = params["initial time"]
    final_time = params["final time"]
    time_step = params["time step"]
    absolute_error = relative_error = 0.0
    time = prev_time = initial_time
    stop = 0
    converged = false
    stop_disp = Vector{Vector{Float64}}(undef, num_domains)
    stop_velo = Vector{Vector{Float64}}(undef, num_domains)
    stop_acce = Vector{Vector{Float64}}(undef, num_domains)
    stop_traction_force = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_disp = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_velo = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_acce = Vector{Vector{Float64}}(undef, num_domains)
    time_hist = Vector{Vector{Float64}}(undef, num_domains)
    disp_hist = Vector{Vector{Vector{Float64}}}(undef, num_domains)
    velo_hist = Vector{Vector{Vector{Float64}}}(undef, num_domains)
    acce_hist = Vector{Vector{Vector{Float64}}}(undef, num_domains)
    traction_force_hist = Vector{Vector{Vector{Float64}}}(undef, num_domains)
    SolidSchwarzController(num_domains, minimum_iterations, maximum_iterations,
        absolute_tolerance, relative_tolerance, absolute_error, relative_error,
        initial_time, final_time, time_step, time, prev_time, stop, converged,
        stop_disp, stop_velo, stop_acce, stop_traction_force, schwarz_disp, schwarz_velo, schwarz_acce,
        time_hist, disp_hist, velo_hist, acce_hist, traction_force_hist)
end

function create_schwarz_controller(params::Dict{Any,Any})
    type = params["subdomains type"]
    if type == "static solid mechanics" || type == "dynamic solid mechanics"
        return SolidSchwarzController(params)
    else
        error("Unknown type of Schwarz controller : ", type)
    end
end

function schwarz(sim::MultiDomainSimulation)
    set_subcycle_times(sim)
    iteration_number = 1
    save_stop_solutions(sim)
    while true
        println("Schwarz iteration=", iteration_number)
        save_schwarz_solutions(sim)
        subcycle(sim)
        iteration_number += 1
        ΔX = update_schwarz_convergence_criterion(sim)
        println("Schwarz criterion |ΔX|=", ΔX)
        if stop_schwarz(sim, iteration_number) == true
            break
        end
        restore_stop_solutions(sim)
    end
end

function save_stop_solutions(sim::MultiDomainSimulation)
    save_stop_solutions(sim.schwarz_controller, sim.subsims)
end

function save_stop_solutions(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    for i ∈ 1:schwarz_controller.num_domains
        schwarz_controller.stop_disp[i] = sims[i].integrator.displacement
        schwarz_controller.stop_velo[i] = sims[i].integrator.velocity
        schwarz_controller.stop_acce[i] = sims[i].integrator.acceleration
        schwarz_controller.stop_traction_force[i] = sims[i].model.boundary_traction_force
    end
end

function restore_stop_solutions(sim::MultiDomainSimulation)
    restore_stop_solutions(sim.schwarz_controller, sim.subsims)
end

function restore_stop_solutions(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    for i ∈ 1:schwarz_controller.num_domains
        sims[i].integrator.displacement = schwarz_controller.stop_disp[i]
        sims[i].integrator.velocity = schwarz_controller.stop_velo[i]
        sims[i].integrator.acceleration = schwarz_controller.stop_acce[i]
        sims[i].model.boundary_traction_force = schwarz_controller.stop_traction_force[i]
        copy_solution_source_targets(sims[i].integrator, sims[i].solver, sims[i].model)
    end
end

function save_schwarz_solutions(sim::MultiDomainSimulation)
    save_schwarz_solutions(sim.schwarz_controller, sim.subsims)
end

function save_schwarz_solutions(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    for i ∈ 1:schwarz_controller.num_domains
        schwarz_controller.schwarz_disp[i] = sims[i].integrator.displacement
        schwarz_controller.schwarz_velo[i] = sims[i].integrator.velocity
        schwarz_controller.schwarz_acce[i] = sims[i].integrator.acceleration
    end
end

function set_subcycle_times(sim::MultiDomainSimulation)
    initial_time = sim.schwarz_controller.prev_time
    final_time = sim.schwarz_controller.time
    for subsim ∈ sim.subsims
        subsim.integrator.initial_time = initial_time
        subsim.integrator.time = initial_time
        subsim.integrator.final_time = final_time
    end
end

function subcycle(sim::MultiDomainSimulation)
    setup_subcycle(sim)
    subsim_index = 1
    for subsim ∈ sim.subsims
        println("subcycle ", subsim.name)
        stop_index = 1
        while true
            advance_time(subsim)
            if stop_evolve(subsim) == true
                break
            end
            subsim.model.time = subsim.integrator.time
            apply_bcs(subsim)
            advance(subsim)
            stop_index += 1
            save_history_snapshot(schwarz_controller, sim.subsims, subsim_index, stop_index)
        end
        subsim_index +=1
    end
end

function setup_subcycle(sim::MultiDomainSimulation)
    resize_histories(sim.schwarz_controller, sim.subsims)
end

function resize_histories(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    for subsim ∈ 1:schwarz_controller.num_domains
        num_steps = round(Int64, sims[subsim].integrator.time_step / schwarz_controller.time_step)
        Δt = schwarz_controller.time_step / num_steps
        num_stops = num_steps + 1
        sims[subsim].integrator.time_step = Δt
        resize!(schwarz_controller.time_hist[subsim], num_stops)
        resize!(schwarz_controller.disp_hist[subsim], num_stops)
        resize!(schwarz_controller.velo_hist[subsim], num_stops)
        resize!(schwarz_controller.acce_hist[subsim], num_stops)
        resize!(schwarz_controller.traction_force_hist[subsim], num_stops)
        for stop ∈ 1:num_stops
            schwarz_controller.time_hist[subsim][stop] = schwarz_controller.prev_time + (stop - 1) * Δt
            schwarz_controller.disp_hist[subsim][stop] = schwarz_controller.stop_disp[subsim]
            schwarz_controller.velo_hist[subsim][stop] = schwarz_controller.stop_velo[subsim]
            schwarz_controller.acce_hist[subsim][stop] = schwarz_controller.stop_acce[subsim]
            schwarz_controller.traction_force_hist[subsim][stop] = schwarz_controller.stop_traction_force[subsim]
        end
    end
end

function save_history_snapshot(schwarz_controller::SchwarzController, sims::Vector{SingleDomainSimulation}, subsim_index::Int64, stop_index::Int64)
    schwarz_controller.disp_hist[subsim_index][stop_index] = sims[subsim_index].integrator.displacement
    schwarz_controller.velo_hist[subsim_index][stop_index] = sims[subsim_index].integrator.velocity
    schwarz_controller.acce_hist[subsim_index][stop_index] = sims[subsim_index].integrator.acceleration
    schwarz_controller.traction_force_hist[subsim_index][stop_index] = sims[subsim_index].model.boundary_traction_force
end

function update_schwarz_convergence_criterion(sim::MultiDomainSimulation)
    return update_schwarz_convergence_criterion(sim.schwarz_controller, sim.subsims)
end

function update_schwarz_convergence_criterion(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    num_sims = length(sims)
    norms_disp = zeros(num_sims)
    norms_diff = zeros(num_sims)
    for i ∈ 1:num_sims
        Δt = schwarz_controller.time_step
        xᵖʳᵉᵛ = schwarz_controller.schwarz_disp[i] + Δt * schwarz_controller.schwarz_velo[i]
        xᶜᵘʳʳ = sims[i].integrator.displacement + Δt * sims[i].integrator.velocity
        norms_disp[i] = norm(xᵖʳᵉᵛ)
        norms_diff[i] = norm(xᶜᵘʳʳ - xᵖʳᵉᵛ)
    end
    norm_disp = norm(norms_disp)
    norm_diff = norm(norms_diff)
    schwarz_controller.absolute_error = norm_diff
    schwarz_controller.relative_error = norm_disp > 0.0 ? norm_diff / norm_disp : norm_diff
    conv_abs = schwarz_controller.absolute_error ≤ schwarz_controller.absolute_tolerance
    conv_rel = schwarz_controller.relative_error ≤ schwarz_controller.relative_tolerance
    schwarz_controller.converged = conv_abs || conv_rel
    return norm_diff
end

function stop_schwarz(sim::MultiDomainSimulation, iteration_number::Int64)
    for subsim ∈ sim.subsims
        if subsim.solver.failed == true
            return true
        end
    end
    if sim.schwarz_controller.absolute_error == 0.0
        return true
    end
    exceeds_minimum_iterations = iteration_number > sim.schwarz_controller.minimum_iterations
    if exceeds_minimum_iterations == false
        return false
    end
    exceeds_maximum_iterations = iteration_number > sim.schwarz_controller.maximum_iterations
    if exceeds_maximum_iterations == true
        return true
    end
    return sim.schwarz_controller.converged
end