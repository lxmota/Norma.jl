@variables t, x, y, z
D = Differential(t)

function SMDirichletBC(input_mesh::PyObject, bc_params::Dict{Any,Any})
    node_set_name = bc_params["node set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = input_mesh.get_node_set_nodes(node_set_id)
    # expression is an arbitrary function of t, x, y, z in the input file
    disp_num = eval(Meta.parse(expression))
    velo_num = expand_derivatives(D(disp_num))
    acce_num = expand_derivatives(D(velo_num))
    SMDirichletBC(node_set_name, offset, node_set_id, node_set_node_indices,
        disp_num, velo_num, acce_num)
end

function SMNeumannBC(input_mesh::PyObject, bc_params::Dict{Any,Any})
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = input_mesh.get_side_set_node_list(side_set_id)
    # expression is an arbitrary function of t, x, y, z in the input file
    traction_num = eval(Meta.parse(expression))
    SMNeumannBC(side_set_name, offset, side_set_id, num_nodes_per_side, side_set_node_indices, traction_num)
end

function SMContactSchwarzBC(coupled_subsim::SingleDomainSimulation, input_mesh::PyObject, bc_params::Dict{Any,Any})
    side_set_name = bc_params["side set"]
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = input_mesh.get_side_set_node_list(side_set_id)
    coupled_block_name = bc_params["source block"]
    coupled_bc_index = 0
    coupled_mesh = coupled_subsim.params["input_mesh"]
    coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
    coupled_side_set_name = bc_params["source side set"]
    coupled_side_set_id = side_set_id_from_name(coupled_side_set_name, coupled_mesh)
    is_dirichlet = false
    SMContactSchwarzBC(side_set_name, side_set_id, num_nodes_per_side, 
        side_set_node_indices, coupled_subsim, coupled_bc_index, coupled_mesh, coupled_block_id, coupled_side_set_id, is_dirichlet)
end

function SMSchwarzDBC(subsim::SingleDomainSimulation, coupled_subsim::SingleDomainSimulation, input_mesh::PyObject, bc_params::Dict{Any,Any})
    node_set_name = bc_params["node set"]
    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = input_mesh.get_node_set_nodes(node_set_id)
    coupled_block_name = bc_params["source block"]
    coupled_mesh = coupled_subsim.params["input_mesh"]
    coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
    element_type = coupled_mesh.elem_type(coupled_block_id)
    coupled_nodes_indices = Vector{Vector{Int64}}(undef, 0)
    interpolation_function_values = Vector{Vector{Float64}}(undef, 0)
    for node_index ∈ node_set_node_indices
        point = subsim.model.reference[:, node_index]
        node_indices, ξ, found = find_in_mesh(point, coupled_subsim.model, coupled_mesh, coupled_block_id)
        if found == false
            error("Could not find point: ", point, " in subdomain: ", coupled_subsim.name)
        end
        N, _, _ = interpolate(element_type, ξ)
        push!(coupled_nodes_indices, node_indices)
        push!(interpolation_function_values, N)
    end
    SMSchwarzDBC(node_set_name, node_set_id, node_set_node_indices,
        coupled_subsim, coupled_mesh, coupled_block_id, coupled_nodes_indices, interpolation_function_values)
end

function apply_bc(model::SolidMechanics, bc::SMDirichletBC)
    for node_index ∈ bc.node_set_node_indices
        values = Dict(t=>model.time, x=>model.reference[1, node_index], y=>model.reference[2, node_index], z=>model.reference[3, node_index])
        disp_sym = substitute(bc.disp_num, values)
        velo_sym = substitute(bc.velo_num, values)
        acce_sym = substitute(bc.acce_num, values)
        disp_val = extract_value(disp_sym)
        velo_val = extract_value(velo_sym)
        acce_val = extract_value(acce_sym)
        dof_index = 3 * (node_index - 1) + bc.offset
        model.current[bc.offset, node_index] = model.reference[bc.offset, node_index] + disp_val
        model.velocity[bc.offset, node_index] = velo_val
        model.acceleration[bc.offset, node_index] = acce_val
        model.free_dofs[dof_index] = false
    end
end

function apply_bc(model::SolidMechanics, bc::SMNeumannBC)
    ss_node_index = 1
    for side ∈ bc.num_nodes_per_side
        side_nodes = bc.side_set_node_indices[ss_node_index:ss_node_index+side-1]
        side_coordinates = model.reference[:, side_nodes]
        nodal_force_component = get_side_set_nodal_forces(side_coordinates, bc.traction_num, model.time)
        ss_node_index += side
        side_node_index = 1
        for node_index ∈ side_nodes
            bc_val = nodal_force_component[side_node_index]
            side_node_index += 1
            dof_index = 3 * (node_index - 1) + bc.offset
            model.boundary_force[dof_index] += bc_val
        end
    end
end

function find_in_mesh(point::Vector{Float64}, model::SolidMechanics, mesh::PyObject, blk_id::Int64)
    element_type = mesh.elem_type(blk_id)
    elem_blk_conn, num_blk_elems, num_elem_nodes = mesh.get_elem_connectivity(blk_id)
    node_indices = Vector{Int64}()
    found = false
    for blk_elem_index ∈ 1:num_blk_elems
        conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
        node_indices = elem_blk_conn[conn_indices]
        elem_ref_pos = model.reference[:, node_indices]
        ξ = map_to_parametric(element_type, elem_ref_pos, point)
        found = is_inside_parametric(element_type, ξ)
        if found == true
            return node_indices, ξ, true
        end
    end
    return node_indices, ξ, false
end

function apply_bc_detail(model::SolidMechanics, bc::SMContactSchwarzBC)
    if bc.is_dirichlet == true
        apply_sm_schwarz_contact_dirichlet(model, bc)
    else
        apply_sm_schwarz_contact_neumann(model, bc)
    end
end

function apply_bc_detail(model::SolidMechanics, bc::SMSchwarzDBC)
    for i ∈ 1:length(bc.node_set_node_indices)
        node_index = bc.node_set_node_indices[i]
        coupled_node_indices = bc.coupled_nodes_indices[i]
        N = bc.interpolation_function_values[i]
        elem_posn = bc.coupled_subsim.model.current[:, coupled_node_indices]
        elem_velo = bc.coupled_subsim.model.velocity[:, coupled_node_indices]
        elem_acce = bc.coupled_subsim.model.acceleration[:, coupled_node_indices]
        point_posn = elem_posn * N
        point_velo = elem_velo * N
        point_acce = elem_acce * N
        model.current[:, node_index] = point_posn
        model.velocity[:, node_index] = point_velo
        model.acceleration[:, node_index] = point_acce
        dof_index = [3 * node_index - 2, 3 * node_index - 1, 3 * node_index]
        model.free_dofs[dof_index] .= false
    end
end

function apply_bc(model::SolidMechanics, bc::SchwarzBoundaryCondition)
    global_sim = bc.coupled_subsim.params["global_simulation"]
    schwarz_controller = global_sim.schwarz_controller
    if schwarz_controller.schwarz_contact == true && schwarz_controller.active_contact == false
        return
    end
    empty_history = length(global_sim.schwarz_controller.time_hist) == 0
    same_step = schwarz_controller.same_step == true
    if empty_history == true
        apply_bc_detail(model, bc)
        return
    end
    # Save solution of coupled simulation
    saved_disp = bc.coupled_subsim.integrator.displacement
    saved_velo = bc.coupled_subsim.integrator.velocity
    saved_acce = bc.coupled_subsim.integrator.acceleration
    saved_∂Ω_f = bc.coupled_subsim.model.internal_force
    time = model.time
    coupled_name = bc.coupled_subsim.name
    coupled_index = global_sim.subsim_name_index_map[coupled_name]
    time_hist = global_sim.schwarz_controller.time_hist[coupled_index]
    disp_hist = global_sim.schwarz_controller.disp_hist[coupled_index]
    velo_hist = global_sim.schwarz_controller.velo_hist[coupled_index]
    acce_hist = global_sim.schwarz_controller.acce_hist[coupled_index]
    ∂Ω_f_hist = global_sim.schwarz_controller.∂Ω_f_hist[coupled_index]
    interp_disp = same_step == true ? disp_hist[end] : interpolate(time_hist, disp_hist, time)
    interp_velo = same_step == true ? velo_hist[end] : interpolate(time_hist, velo_hist, time)
    interp_acce = same_step == true ? acce_hist[end] : interpolate(time_hist, acce_hist, time)
    interp_∂Ω_f = same_step == true ? ∂Ω_f_hist[end] : interpolate(time_hist, ∂Ω_f_hist, time)
    bc.coupled_subsim.integrator.displacement = interp_disp
    bc.coupled_subsim.integrator.velocity = interp_velo
    bc.coupled_subsim.integrator.acceleration = interp_acce
    bc.coupled_subsim.model.internal_force = interp_∂Ω_f
    copy_solution_source_targets(bc.coupled_subsim.integrator, bc.coupled_subsim.solver, bc.coupled_subsim.model)
    apply_bc_detail(model, bc)
    bc.coupled_subsim.integrator.displacement = saved_disp
    bc.coupled_subsim.integrator.velocity = saved_velo
    bc.coupled_subsim.integrator.acceleration = saved_acce
    bc.coupled_subsim.model.internal_force = saved_∂Ω_f
    copy_solution_source_targets(bc.coupled_subsim.integrator, bc.coupled_subsim.solver, bc.coupled_subsim.model)
end

function transfer_normal_component(source::Vector{Float64}, target::Vector{Float64}, normal::Vector{Float64})
    normal_projection = normal * normal'
    tangent_projection = I(length(normal)) - normal_projection
    return tangent_projection * target + normal_projection * source
end

function apply_sm_schwarz_contact_dirichlet(model::SolidMechanics, bc::SMContactSchwarzBC)
    ss_node_index = 1
    for side ∈ bc.num_nodes_per_side
        side_nodes = bc.side_set_node_indices[ss_node_index:ss_node_index+side-1]
        ss_node_index += side
        for node_index ∈ side_nodes
            point = model.current[:, node_index]
            point_new, ξ, closest_face_node_coords, closest_face_node_indices, closest_normal, found = find_and_project(point, bc.coupled_mesh, bc.coupled_side_set_id, bc.coupled_subsim.model)
            if found == false
                is_inside, points_inside, int_points_coords, all_dst_face_node_coords, all_dst_face_node_indices = search_integration_points(side_nodes, model, bc)
                space_dim, num_int_points = size(int_points_coords)
                parametric_dim = space_dim - 1
                if is_inside == true
                    minimum_distance = Inf
                    for int_point ∈ 1:num_int_points
                        int_point_coords = int_points_coords[:, int_point]
                        diff = point - int_point_coords
                        distance = norm(diff)
                        if distance < minimum_distance
                            if points_inside[int_point] == true
                                minimum_distance = distance
                                closest_face_node_coords = all_dst_face_node_coords[int_point]
                                closest_face_node_indices = all_dst_face_node_indices[int_point]
                            end
                        end
                    end
                    point_new, ξ, _, closest_normal = closest_point_projection(parametric_dim, closest_face_node_coords, point)
                end
            end
            if found == true || is_inside == true
                model.current[:, node_index] = point_new
                element_type = get_element_type(2, side)
                N, _, _ = interpolate(element_type, ξ)
                source_velo = bc.coupled_subsim.model.velocity[:, closest_face_node_indices] * N
                source_acce = bc.coupled_subsim.model.acceleration[:, closest_face_node_indices] * N
                model.velocity[:, node_index] = transfer_normal_component(source_velo, model.velocity[:, node_index], closest_normal)
                model.acceleration[:, node_index] = transfer_normal_component(source_acce, model.acceleration[:, node_index], closest_normal)
                dof_index = [3 * node_index - 2]
                model.free_dofs[dof_index] .= false
            else    
                continue
            end
        end
    end
end

function apply_sm_schwarz_contact_neumann(model::SolidMechanics, bc::SMContactSchwarzBC)
    schwarz_tractions, normals = get_dst_traction(model, bc)
    local_to_global_map = get_side_set_local_to_global_map(model.mesh, bc.side_set_id)
    num_local_nodes = length(local_to_global_map)
    for local_node ∈ 1:num_local_nodes
        global_node = local_to_global_map[local_node]
        node_tractions = schwarz_tractions[3*local_node-2:3*local_node]
        normal = normals[:, local_node]
        model.boundary_force[3*global_node-2:3*global_node] += transfer_normal_component(node_tractions, model.boundary_force[3*global_node-2:3*global_node], normal)
    end
end

function reduce_traction(mesh::PyObject, side_set_id::Integer, global_traction::Vector{Float64})
    local_to_global_map = get_side_set_local_to_global_map(mesh, side_set_id)
    num_local_nodes = length(local_to_global_map)
    local_traction = zeros(3*num_local_nodes)
    for local_node ∈ 1:num_local_nodes
        global_node = local_to_global_map[local_node]
        local_traction[3*local_node-2:3*local_node] = global_traction[3*global_node-2:3*global_node]
    end
    return local_traction
end

function get_dst_traction(dst_model::SolidMechanics, bc::SMContactSchwarzBC)
    src_mesh = bc.coupled_subsim.model.mesh
    src_side_set_id = bc.coupled_side_set_id
    src_global_traction = -bc.coupled_subsim.model.internal_force
    src_model = bc.coupled_subsim.model
    dst_mesh = dst_model.mesh
    dst_side_set_id = bc.side_set_id
    square_projection_matrix = get_square_projection_matrix(src_mesh, src_model, src_side_set_id)
    rectangular_projection_matrix, normals = get_rectangular_projection_matrix(dst_mesh, dst_model, dst_side_set_id, src_mesh, src_model, src_side_set_id)
    src_local_traction = reduce_traction(src_mesh, src_side_set_id, src_global_traction)
    src_traction_x = src_local_traction[1:3:end]
    src_traction_y = src_local_traction[2:3:end]
    src_traction_z = src_local_traction[3:3:end]
    projection_operator = rectangular_projection_matrix * inv(square_projection_matrix)
    dst_traction_x = projection_operator * src_traction_x
    dst_traction_y = projection_operator * src_traction_y
    dst_traction_z = projection_operator * src_traction_z
    dst_traction = zeros(3*length(dst_traction_x))
    dst_traction[1:3:end] = dst_traction_x
    dst_traction[2:3:end] = dst_traction_y
    dst_traction[3:3:end] = dst_traction_z
    return dst_traction, normals
end    

function node_set_id_from_name(node_set_name::String, mesh::PyObject)
    node_set_names = mesh.get_node_set_names()
    num_names = length(node_set_names)
    node_set_index = 0
    for index ∈ 1:num_names
        if (node_set_name == node_set_names[index])
            node_set_index = index
            break
        end
    end
    if (node_set_index == 0)
        error("node set ", node_set_name, " cannot be found in mesh")
    end
    node_set_ids = mesh.get_node_set_ids()
    node_set_id = node_set_ids[node_set_index]
    return Int64(node_set_id)
end

function side_set_id_from_name(side_set_name::String, mesh::PyObject)
    side_set_names = mesh.get_side_set_names()
    num_names = length(side_set_names)
    side_set_index = 0
    for index ∈ 1:num_names
        if (side_set_name == side_set_names[index])
            side_set_index = index
            break
        end
    end
    if (side_set_index == 0)
        error("side set ", side_set_name, " cannot be found in mesh")
    end
    side_set_ids = mesh.get_side_set_ids()
    side_set_id = side_set_ids[side_set_index]
    return Int64(side_set_id)
end

function block_id_from_name(block_name::String, mesh::PyObject)
    block_names = mesh.get_elem_blk_names()
    num_names = length(block_names)
    block_index = 0
    for index ∈ 1:num_names
        if (block_name == block_names[index])
            block_index = index
            break
        end
    end
    if (block_index == 0)
        error("block ", block_name, " cannot be found in mesh")
    end
    block_ids = mesh.get_elem_blk_ids()
    block_id = block_ids[block_index]
    return Int64(block_id)
end

function component_offset_from_string(name::String)
    offset = 0
    if name == "x"
        offset = 1
    elseif name == "y"
        offset = 2
    elseif name == "z"
        offset = 3
    else
        error("invalid component name ", name)
    end
    return offset
end

function extract_value(value::Real)
    return value
end

function extract_value(symbol::Num)
    return symbol.val
end

function create_bcs(params::Dict{Any,Any})
    boundary_conditions = Vector{BoundaryCondition}()
    if haskey(params, "boundary conditions") == false
        return boundary_conditions
    end
    input_mesh = params["input_mesh"]
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) ∈ bc_params
        for bc_setting_params ∈ bc_type_params
            if bc_type == "Dirichlet"
                boundary_condition = SMDirichletBC(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Neumann"
                boundary_condition = SMNeumannBC(input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)
            elseif bc_type == "Schwarz contact"
                sim = params["global_simulation"]
                sim.schwarz_controller.schwarz_contact = true
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                boundary_condition = SMContactSchwarzBC(coupled_subsim, input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)                
            elseif bc_type == "Schwarz Dirichlet"
                sim = params["global_simulation"]
                subsim_name = params["name"]
                subdomain_index = sim.subsim_name_index_map[subsim_name]
                subsim = sim.subsims[subdomain_index]
                coupled_subsim_name = bc_setting_params["source"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_subsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                boundary_condition = SMSchwarzDBC(subsim, coupled_subsim, input_mesh, bc_setting_params)
                push!(boundary_conditions, boundary_condition)                
            else
                error("Unknown boundary condition type : ", bc_type)
            end
        end
    end
    return boundary_conditions
end

function apply_bcs(model::SolidMechanics)
    _, num_nodes = size(model.reference)
    model.boundary_force = zeros(3*num_nodes)
    model.free_dofs = trues(3 * num_nodes)
    for boundary_condition ∈ model.boundary_conditions
        apply_bc(model, boundary_condition)
    end
end

function apply_ics(params::Dict{Any,Any}, model::SolidMechanics)
    if haskey(params, "initial conditions") == false
        return
    end
    input_mesh = params["input_mesh"]
    ic_params = params["initial conditions"]
    for (ic_type, ic_type_params) ∈ ic_params
        for ic ∈ ic_type_params
            node_set_name = ic["node set"]
            expr_str = ic["function"]
            component = ic["component"]
            offset = component_offset_from_string(component)
            node_set_id = node_set_id_from_name(node_set_name, input_mesh)
            node_set_node_indices = input_mesh.get_node_set_nodes(node_set_id)
            # expr_str is an arbitrary function of x, y, z in the input file
            ic_expr = Meta.parse(expr_str)
            ic_eval = eval(ic_expr)
            for node_index ∈ node_set_node_indices
                values = Dict(x=>model.reference[1, node_index], y=>model.reference[2, node_index], z=>model.reference[3, node_index])
                ic_sym = substitute(ic_eval, values)
                ic_val = extract_value(ic_sym)
                if ic_type == "displacement"
                    model.current[offset, node_index] = model.reference[offset, node_index] + ic_val
                elseif ic_type == "velocity"
                    model.velocity[offset, node_index] = ic_val
                end
            end
        end
    end
end

function pair_schwarz_bcs(sim::MultiDomainSimulation)
    for subsim ∈ sim.subsims
        model = subsim.model
        name = subsim.name
        bcs = model.boundary_conditions
        for bc ∈ bcs
            pair_bc(name, bc)
        end
    end
end

function pair_bc(_::String, _::RegularBoundaryCondition)
end

function pair_bc(_::String, _::SchwarzBoundaryCondition)
end

function pair_bc(name::String, bc::ContactSchwarzBoundaryCondition)
    coupled_model = bc.coupled_subsim.model
    coupled_bcs = coupled_model.boundary_conditions
    for coupled_bc ∈ coupled_bcs
        if is_coupled_to_current(name, coupled_bc) == true
            coupled_bc.is_dirichlet = !bc.is_dirichlet
        end
    end
end

function is_coupled_to_current(_::String, _::RegularBoundaryCondition)
    return false
end

function is_coupled_to_current(_::String, _::SchwarzBoundaryCondition)
    return false
end

function is_coupled_to_current(name::String, coupled_bc::ContactSchwarzBoundaryCondition)
    return name == coupled_bc.coupled_subsim.name
end