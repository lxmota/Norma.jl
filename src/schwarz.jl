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
    schwarz_disp = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_velo = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_acce = Vector{Vector{Float64}}(undef, num_domains)
    time_hist = Vector{Float64}(undef, num_domains)
    disp_hist = Vector{Vector{Vector{Float64}}}(undef, num_domains)
    velo_hist = Vector{Vector{Vector{Float64}}}(undef, num_domains)
    acce_hist = Vector{Vector{Vector{Float64}}}(undef, num_domains)
    SolidSchwarzController(num_domains, minimum_iterations, maximum_iterations,
        absolute_tolerance, relative_tolerance, absolute_error, relative_error,
        initial_time, final_time, time_step, time, prev_time, stop, converged,
        stop_disp, stop_velo, stop_acce, schwarz_disp, schwarz_velo, schwarz_acce,
        time_hist, disp_hist, velo_hist, acce_hist)
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
    for subsim ∈ sim.subsims
        println("subcycle ", subsim.name)
        while true
            advance_time(subsim)
            if stop_evolve(subsim) == true
                break
            end
            subsim.model.time = subsim.integrator.time
            apply_bcs(subsim)
            advance(subsim)
        end
    end
end

function setup_subcycle(sim::MultiDomainSimulation)
    resize_histories(sim.schwarz_controller, sim.subsims)
end

function resize_histories(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    for i ∈ 1:schwarz_controller.num_domains
        num_steps = round(Int64, sims[i].integrator.time_step / schwarz_controller.time_step)
        sims[i].integrator.time_step = schwarz_controller.time_step / num_steps
        resize!(schwarz_controller.disp_hist[i], num_steps)
        resize!(schwarz_controller.velo_hist[i], num_steps)
        resize!(schwarz_controller.acce_hist[i], num_steps)
        for j ∈ 1:num_steps
            schwarz_controller.disp_hist[i][j] = schwarz_controller.stop_disp[i]
            schwarz_controller.velo_hist[i][j] = schwarz_controller.stop_velo[i]
            schwarz_controller.acce_hist[i][j] = schwarz_controller.stop_acce[i]
        end
    end
end

function update_schwarz_convergence_criterion(sim::MultiDomainSimulation)
    return is_schwarz_converged(sim.schwarz_controller, sim.subsims)
end

function is_schwarz_converged(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
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