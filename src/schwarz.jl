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
    same_step = get(params, "same time step for domains", false)
    stop = 0
    converged = false
    stop_disp = Vector{Vector{Float64}}(undef, num_domains)
    stop_velo = Vector{Vector{Float64}}(undef, num_domains)
    stop_acce = Vector{Vector{Float64}}(undef, num_domains)
    stop_traction_force = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_disp = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_velo = Vector{Vector{Float64}}(undef, num_domains)
    schwarz_acce = Vector{Vector{Float64}}(undef, num_domains)
    time_hist = Vector{Vector{Float64}}()
    disp_hist = Vector{Vector{Vector{Float64}}}()
    velo_hist = Vector{Vector{Vector{Float64}}}()
    acce_hist = Vector{Vector{Vector{Float64}}}()
    traction_force_hist = Vector{Vector{Vector{Float64}}}()
    schwarz_contact = false
    active_contact = false
    SolidSchwarzController(num_domains, minimum_iterations, maximum_iterations,
        absolute_tolerance, relative_tolerance, absolute_error, relative_error,
        initial_time, final_time, time_step, time, prev_time, same_step, stop, converged,
        stop_disp, stop_velo, stop_acce, stop_traction_force, schwarz_disp, schwarz_velo, schwarz_acce,
        time_hist, disp_hist, velo_hist, acce_hist, traction_force_hist, schwarz_contact, active_contact)
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
    iteration_number = 1
    save_stop_solutions(sim)
    save_schwarz_solutions(sim)
    resize_histories(sim)
    while true
        println("Schwarz iteration=", iteration_number)
        synchronize(sim)
        set_subcycle_times(sim)
        subcycle(sim)
        iteration_number += 1
        ΔX = update_schwarz_convergence_criterion(sim)
        println("Schwarz criterion |ΔX|=", ΔX)
        if stop_schwarz(sim, iteration_number) == true
            break
        end
        write_step(sim)
        save_schwarz_solutions(sim)
        restore_stop_solutions(sim)
    end
end

function save_stop_solutions(sim::MultiDomainSimulation)
    save_stop_solutions(sim.schwarz_controller, sim.subsims)
end

function save_stop_solutions(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    for i ∈ 1:schwarz_controller.num_domains
        schwarz_controller.stop_disp[i] = deepcopy(sims[i].integrator.displacement)
        schwarz_controller.stop_velo[i] = deepcopy(sims[i].integrator.velocity)
        schwarz_controller.stop_acce[i] = deepcopy(sims[i].integrator.acceleration)
        schwarz_controller.stop_traction_force[i] = deepcopy(sims[i].model.internal_force)
    end
end

function restore_stop_solutions(sim::MultiDomainSimulation)
    restore_stop_solutions(sim.schwarz_controller, sim.subsims)
end

function restore_stop_solutions(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    for i ∈ 1:schwarz_controller.num_domains
        sims[i].integrator.displacement = deepcopy(schwarz_controller.stop_disp[i])
        sims[i].integrator.velocity = deepcopy(schwarz_controller.stop_velo[i])
        sims[i].integrator.acceleration = deepcopy(schwarz_controller.stop_acce[i])
        sims[i].model.internal_force = deepcopy(schwarz_controller.stop_traction_force[i])
        copy_solution_source_targets(sims[i].integrator, sims[i].solver, sims[i].model)
    end
end

function save_schwarz_solutions(sim::MultiDomainSimulation)
    save_schwarz_solutions(sim.schwarz_controller, sim.subsims)
end

function save_schwarz_solutions(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    for i ∈ 1:schwarz_controller.num_domains
        schwarz_controller.schwarz_disp[i] = deepcopy(sims[i].integrator.displacement)
        schwarz_controller.schwarz_velo[i] = deepcopy(sims[i].integrator.velocity)
        schwarz_controller.schwarz_acce[i] = deepcopy(sims[i].integrator.acceleration)
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
            save_history_snapshot(sim.schwarz_controller, sim.subsims, subsim_index, stop_index)
        end
        subsim_index +=1
    end
end

function resize_histories(sim::MultiDomainSimulation)
    resize_histories(sim.schwarz_controller, sim.subsims)
end

function resize_histories(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    num_domains = schwarz_controller.num_domains
    resize!(schwarz_controller.time_hist, num_domains)
    resize!(schwarz_controller.disp_hist, num_domains)
    resize!(schwarz_controller.velo_hist, num_domains)
    resize!(schwarz_controller.acce_hist, num_domains)
    resize!(schwarz_controller.traction_force_hist, num_domains)
    for subsim ∈ 1:num_domains
        num_steps = round(Int64, schwarz_controller.time_step / sims[subsim].integrator.time_step)
        Δt = schwarz_controller.time_step / num_steps
        num_stops = num_steps + 1
        sims[subsim].integrator.time_step = Δt
        schwarz_controller.time_hist[subsim] = Vector{Float64}(undef, num_stops)
        schwarz_controller.disp_hist[subsim] = Vector{Vector{Float64}}(undef, num_stops)
        schwarz_controller.velo_hist[subsim] = Vector{Vector{Float64}}(undef, num_stops)
        schwarz_controller.acce_hist[subsim] = Vector{Vector{Float64}}(undef, num_stops)
        schwarz_controller.traction_force_hist[subsim] = Vector{Vector{Float64}}(undef, num_stops)
        for stop ∈ 1:num_stops
            schwarz_controller.time_hist[subsim][stop] = schwarz_controller.prev_time + (stop - 1) * Δt
            schwarz_controller.disp_hist[subsim][stop] = deepcopy(schwarz_controller.stop_disp[subsim])
            schwarz_controller.velo_hist[subsim][stop] = deepcopy(schwarz_controller.stop_velo[subsim])
            schwarz_controller.acce_hist[subsim][stop] = deepcopy(schwarz_controller.stop_acce[subsim])
            schwarz_controller.traction_force_hist[subsim][stop] = deepcopy(schwarz_controller.stop_traction_force[subsim])
        end
    end
end

function save_history_snapshot(schwarz_controller::SchwarzController, sims::Vector{SingleDomainSimulation}, subsim_index::Int64, stop_index::Int64)
    schwarz_controller.disp_hist[subsim_index][stop_index] = deepcopy(sims[subsim_index].integrator.displacement)
    schwarz_controller.velo_hist[subsim_index][stop_index] = deepcopy(sims[subsim_index].integrator.velocity)
    schwarz_controller.acce_hist[subsim_index][stop_index] = deepcopy(sims[subsim_index].integrator.acceleration)
    schwarz_controller.traction_force_hist[subsim_index][stop_index] = deepcopy(sims[subsim_index].model.internal_force)
end

function update_schwarz_convergence_criterion(sim::MultiDomainSimulation)
    return update_schwarz_convergence_criterion(sim.schwarz_controller, sim.subsims)
end

function update_schwarz_convergence_criterion(schwarz_controller::SolidSchwarzController, sims::Vector{SingleDomainSimulation})
    num_domains = schwarz_controller.num_domains
    norms_disp = zeros(num_domains)
    norms_diff = zeros(num_domains)
    for i ∈ 1:num_domains
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

function detect_contact(_::SingleDomainSimulation)
end

function detect_contact(sim::MultiDomainSimulation)
    if sim.schwarz_controller.schwarz_contact == false
        return
    end
    num_domains = sim.schwarz_controller.num_domains
    #first condition: contact at the previous stop
    contact_prev = sim.schwarz_controller.active_contact
    contact_domain = zeros(Bool, num_domains)
    for i ∈ 1:sim.schwarz_controller.num_domains
        subsim = sim.subsims[i]
        subsim_bc_params = subsim.params["boundary conditions"]
        mesh = subsim.params["input_mesh"]
        Schwarz_bcs = subsim_bc_params["Schwarz contact"]
        for bc_setting_params ∈ Schwarz_bcs
            side_set_name = bc_setting_params["side set"]
            side_set_id = side_set_id_from_name(side_set_name, mesh)
            num_nodes_per_side, side_set_node_indices = mesh.get_side_set_node_list(side_set_id)
            coupled_subsim_name = bc_setting_params["source"]
            coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
            coupled_subsim = sim.subsims[coupled_subdomain_index]
            coupled_mesh = coupled_subsim.params["input_mesh"]
            coupled_side_set_name = bc_setting_params["source side set"]
            coupled_side_set_id = side_set_id_from_name(coupled_side_set_name, coupled_mesh)
            #second condition: overlap
            global_to_local_map, _, _ = get_side_set_global_to_local_map(mesh, side_set_id)
            overlap_nodes = zeros(Bool,length(global_to_local_map))
            ss_node_index = 1
            for side ∈ num_nodes_per_side
                side_nodes = side_set_node_indices[ss_node_index:ss_node_index+side-1]
                for node_index ∈ side_nodes
                    point = subsim.model.current[:, node_index]
                    _, _, _, _, _, found = find_and_project(point, coupled_mesh, coupled_side_set_id, coupled_subsim.model)
                    node = get.(Ref(global_to_local_map), node_index, 0)
                    overlap_nodes[node] = found
                end
                ss_node_index += side
            end
            #true if at least one node interpenetrates the coupled domain
            overlap = any(overlap_nodes)
            #third condition: compression 
            coupled_global_traction = coupled_subsim.model.internal_force
            coupled_local_traction = reduce_traction(coupled_mesh, coupled_side_set_id, coupled_global_traction)
            compress_nodes = coupled_local_traction[isless.(coupled_local_traction, 0)]
            #true if at least one node interpenetrates the coupled domain
            compress = length(compress_nodes) > 1
            persist = compress && contact_prev
            contact_domain[i] = overlap || persist
        end
    end
    sim.schwarz_controller.active_contact = any(contact_domain)
    return sim.schwarz_controller.active_contact
end