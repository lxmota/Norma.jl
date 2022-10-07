@enum DOF free Dirichlet Neumann Schwarz

function node_set_id_from_name(node_set_name::String, mesh_struct::PyObject)
    node_set_names = mesh_struct.get_node_set_names()
    num_names = length(node_set_names)
    node_set_index = 0
    for index ∈ 1 : num_names
        if (node_set_name == node_set_names[index])
            node_set_index = index
            break
        end
    end
    if (node_set_index == 0)
        error("node set ", node_set_name, " cannot be found in mesh")
    end
    node_set_ids = mesh_struct.get_node_set_ids()
    node_set_id = node_set_ids[node_set_index]
    return node_set_id
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

function apply_bcs(model::SolidMechanics)
    params = model.params
    reference = model.reference
    current = model.current
    mesh_struct = params["mesh_struct"]
    global t = model.time
    xc, yc, zc = mesh_struct.get_coords()
    num_nodes = length(xc)
    nodal_dofs = [free::DOF for _ ∈ 1 : 3 * num_nodes]
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) ∈ bc_params
        for bc ∈ bc_type_params
            if bc_type == "Dirichlet"
                node_set_name = bc["node set"]
                function_str = bc["function"]
                component = bc["component"]
                offset = component_offset_from_string(component)
                node_set_id = node_set_id_from_name(node_set_name, mesh_struct)
                node_set_node_indices = mesh_struct.get_node_set_nodes(node_set_id)
                for node_index ∈ node_set_node_indices
                    global x = xc[node_index]
                    global y = yc[node_index]
                    global z = zc[node_index]
                    # function_str is an arbitrary function of t, x, y, z in the input file
                    bc_expr = Meta.parse(function_str)
                    bc_val = eval(bc_expr)
                    dof_index = 3 * (node_index - 1) + offset
                    current[dof_index] = reference[dof_index] + bc_val
                    nodal_dofs[dof_index] = Dirichlet::DOF
                end
            elseif bc_type == "Schwarz"
            elseif bc_type == "Neumann"
            end
        end
    end
end

function apply_bcs(model::HeatConduction)
    params = model.params
    temperature = model.temperature
    mesh_struct = params["mesh_struct"]
    global t = model.time
    xc, yc, zc = mesh_struct.get_coords()
    num_nodes = length(xc)
    nodal_dofs = [free::DOF for _ ∈ 1 : num_nodes]
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) ∈ bc_params
        for bc ∈ bc_type_params
            if bc_type == "Dirichlet"
                node_set_name = bc["node set"]
                function_str = bc["function"]
                node_set_id = node_set_id_from_name(node_set_name, mesh_struct)
                node_set_node_indices = mesh_struct.get_node_set_nodes(node_set_id)
                for node_index ∈ node_set_node_indices
                    global x = xc[node_index]
                    global y = yc[node_index]
                    global z = zc[node_index]
                    # function_str is an arbitrary function of t, x, y, z in the input file
                    bc_expr = Meta.parse(function_str)
                    bc_val = eval(bc_expr)
                    temperature[node_index] = bc_val
                    nodal_dofs[node_index] = Dirichlet::DOF
                end
            elseif bc_type == "Schwarz"
            elseif bc_type == "Neumann"
            end
        end
    end
end