include("bcs_def.jl")

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
    return node_set_id
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
    return side_set_id
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
    return block_id
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

using Symbolics
@variables t, x, y, z
D = Differential(t)

function extract_value(value::Real)
    return value
end

function extract_value(symbol::Num)
    return symbol.val
end

function apply_bcs(params::Dict{Any,Any}, model::SolidMechanics)
    if haskey(params, "boundary conditions") == false
        return
    end
    input_mesh = params["input_mesh"]
    _, num_nodes = size(model.reference)
    model.boundary_tractions_force = zeros(3*num_nodes)
    model.free_dofs = trues(3 * num_nodes)
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) ∈ bc_params
        for bc ∈ bc_type_params
            if bc_type == "Dirichlet"
                node_set_name = bc["node set"]
                expr_str = bc["function"]
                component = bc["component"]
                offset = component_offset_from_string(component)
                node_set_id = node_set_id_from_name(node_set_name, input_mesh)
                node_set_node_indices = input_mesh.get_node_set_nodes(node_set_id)
                # expr_str is an arbitrary function of t, x, y, z in the input file
                bc_expr = Meta.parse(expr_str)
                disp_eval = eval(bc_expr)
                velo_eval = expand_derivatives(D(disp_eval))
                acce_eval = expand_derivatives(D(velo_eval))
                for node_index ∈ node_set_node_indices
                    values = Dict(t=>model.time, x=>model.reference[1, node_index], y=>model.reference[2, node_index], z=>model.reference[3, node_index])
                    disp_sym = substitute(disp_eval, values)
                    velo_sym = substitute(velo_eval, values)
                    acce_sym = substitute(acce_eval, values)
                    disp_val = extract_value(disp_sym)
                    velo_val = extract_value(velo_sym)
                    acce_val = extract_value(acce_sym)
                    dof_index = 3 * (node_index - 1) + offset
                    model.current[offset, node_index] = model.reference[offset, node_index] + disp_val
                    model.velocity[offset, node_index] = velo_val
                    model.acceleration[offset, node_index] = acce_val
                    model.free_dofs[dof_index] = false
                end
            elseif bc_type == "Neumann"
                side_set_name = bc["side set"]
                expr_str = bc["function"]
                component = bc["component"]
                offset = component_offset_from_string(component)
                side_set_id = side_set_id_from_name(side_set_name, input_mesh)
                ss_num_nodes_per_side, ss_nodes = input_mesh.get_side_set_node_list(side_set_id)
                # expr_str is an arbitrary function of t, x, y, z in the input file
                bc_expr = Meta.parse(expr_str)
                ss_node_index = 1
                for side ∈ ss_num_nodes_per_side
                    side_nodes = ss_nodes[ss_node_index:ss_node_index+side-1]
                    side_coordinates = model.reference[:, side_nodes]
                    nodal_force_component = get_side_set_nodal_forces(side_coordinates, bc_expr, model.time)
                    ss_node_index += side
                    side_node_index = 1
                    for node_index ∈ side_nodes
                        bc_val = nodal_force_component[side_node_index]
                        side_node_index += 1
                        dof_index = 3 * (node_index - 1) + offset
                        model.boundary_tractions_force[dof_index] += bc_val
                    end
                end
            elseif bc_type == "Schwarz contact Dirichlet"
                side_set_name = bc["side set"]
                component = bc["component"]
                side_set_id = side_set_id_from_name(side_set_name, input_mesh)
                ss_num_nodes_per_side, ss_nodes = input_mesh.get_side_set_node_list(side_set_id)
                #getting the coupled mesh
                coupled_cubsim_name = bc["source"]
                sim = params["global_simulation"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_cubsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                coupled_mesh = coupled_subsim.params["input_mesh"]
                coupled_block_name = bc["source block"]
                coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
                coupled_side_set = bc["source side set"]
                #loop over nodes on the contact side set to get Schwarz displacements
                for side ∈ ss_num_nodes_per_side
                    side_nodes = ss_nodes[ss_node_index:ss_node_index+side-1]
                    for node_index ∈ side_nodes
                        point = current[:, node_index]
                        find_element_for_transfer(point, coupled_mesh, coupled_block_id, coupled_side_set, model)
                    end
                end    
            elseif bc_type == "Schwarz contact Neumann"
                side_set_name = bc["side set"]
                component = bc["component"]
                side_set_id = side_set_id_from_name(side_set_name, input_mesh)
                ss_num_nodes_per_side, ss_nodes = input_mesh.get_side_set_node_list(side_set_id)
                #getting the coupled mesh
                coupled_cubsim_name = bc["source"]
                sim = params["global_simulation"]
                coupled_subdomain_index = sim.subsim_name_index_map[coupled_cubsim_name]
                coupled_subsim = sim.subsims[coupled_subdomain_index]
                coupled_mesh = coupled_subsim.params["input_mesh"]
                coupled_block_name = bc["source block"]
                coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
                coupled_side_set = bc["source side set"]
                #loop over nodes on the contact side set to get Schwarz traction
                for side ∈ ss_num_nodes_per_side
                    side_nodes = ss_nodes[ss_node_index:ss_node_index+side-1]
                    for node_index ∈ side_nodes
                        point = current[:, node_index]
                    end
                end    
            end
        end
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