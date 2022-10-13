using SparseArrays

include("constitutive.jl")

abstract type Model end
abstract type SolidMechanics <: Model end
abstract type HeatConduction <: Model end

@enum DOF free Dirichlet Neumann Schwarz

mutable struct StaticSolid <: SolidMechanics
    params::Dict{Any, Any}
    materials::Vector{Solid}
    reference::Matrix{MTScalar}
    current::Matrix{MTScalar}
    nodal_dofs::Vector{DOF}
    time::MTScalar
end

function StaticSolid(params::Dict{Any, Any})
    mesh_struct = params["mesh_struct"]
    x, y, z = mesh_struct.get_coords()
    num_nodes = length(x)
    reference = Matrix{MTScalar}(undef, 3, num_nodes)
    current = Matrix{MTScalar}(undef, 3, num_nodes)
    for node ∈ 1 : num_nodes
        reference[:, node] = [x[node], y[node], z[node]]
        current[:, node] = [x[node], y[node], z[node]]
    end
    materials_file = params["material"]
    material_params = YAML.load_file(materials_file)
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    if (num_blks_params ≠ num_blks)
        error("number of blocks in mesh ", params["mesh"], " (", num_blks,
        ") must be equal to number of blocks in materials file ", params["material"],
        " (", num_blks_params, ")")
    end
    elem_blk_names = mesh_struct.get_elem_blk_names()
    materials = Vector{Solid}(undef, 0)
    for elem_blk_name ∈ elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        push!(materials, material_model)
    end
    time = 0.0
    nodal_dofs = [free::DOF for _ ∈ 1 : 3 * num_nodes]
    StaticSolid(params, materials, reference, current, nodal_dofs, time)
end

mutable struct DynamicSolid <: SolidMechanics
    params::Dict{Any, Any}
    materials::Vector{Solid}
    reference::Matrix{MTScalar}
    current::Matrix{MTScalar}
    velocity::Matrix{MTScalar}
    acceleration::Matrix{MTScalar}
    nodal_dofs::Vector{DOF}
    time::MTScalar
end

function DynamicSolid(params::Dict{Any, Any})
    mesh_struct = params["mesh_struct"]
    x, y, z = mesh_struct.get_coords()
    num_nodes = length(x)
    reference = Matrix{MTScalar}(undef, 3, num_nodes)
    current = Matrix{MTScalar}(undef, 3, num_nodes)
    for node ∈ 1 : num_nodes
        reference[:, node] = [x[node], y[node], z[node]]
        current[:, node] = [x[node], y[node], z[node]]
        velocity[:, node] = [0.0, 0.0, 0.0]
        acceleration[:, node] = [0.0, 0.0, 0.0]
    end
    materials_file = params["material"]
    material_params = YAML.load_file(materials_file)
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    if (num_blks_params ≠ num_blks)
        error("number of blocks in mesh ", params["mesh"], " (", num_blks,
        ") must be equal to number of blocks in materials file ", params["material"],
        " (", num_blks_params, ")")
    end
    elem_blk_names = mesh_struct.get_elem_blk_names()
    materials = Vector{Solid}(undef, 0)
    for elem_blk_name ∈ elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        push!(materials, material_model)
    end
    time = 0.0
    nodal_dofs = [free::DOF for _ ∈ 1 : 3 * num_nodes]
    DynamicSolid(params, materials, reference, current, velocity, acceleration, nodal_dofs, time)
end

mutable struct StaticHeat <: HeatConduction
    params::Dict{Any, Any}
    materials::Vector{Thermal}
    reference::Matrix{MTScalar}
    temperature::Vector{MTScalar}
    nodal_dofs::Vector{DOF}
    time::MTScalar
end

function StaticHeat(params::Dict{Any, Any})
    mesh_struct = params["mesh_struct"]
    x, y, z = mesh_struct.get_coords()
    num_nodes = length(x)
    reference = Matrix{MTScalar}(undef, 3, num_nodes)
    temperature = Vector{MTScalar}(undef, num_nodes)
    for node ∈ 1 : num_nodes
        reference[:, node] = [x[node], y[node], z[node]]
        temperature[node] = 0.0
    end
    materials_file = params["material"]
    material_params = YAML.load_file(materials_file)
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    if (num_blks_params ≠ num_blks)
        error("number of blocks in mesh ", params["mesh"], " (", num_blks,
        ") must be equal to number of blocks in materials file ", params["material"],
        " (", num_blks_params, ")")
    end
    elem_blk_names = mesh_struct.get_elem_blk_names()
    materials = Vector{Thermal}(undef, 0)
    for elem_blk_name ∈ elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        push!(materials, material_model)
    end
    time = 0.0
    nodal_dofs = [free::DOF for _ ∈ 1 : num_nodes]
    StaticHeat(params, materials, reference, temperature, nodal_dofs, time)
end

mutable struct DynamicHeat <: HeatConduction
    params::Dict{Any, Any}
    materials::Vector{Vector}
    reference::Matrix{MTScalar}
    temperature::Vector{MTScalar}
    rate::Vector{MTScalar}
    nodal_dofs::Vector{DOF}
    time::MTScalar
end

function DynamicHeat(params::Dict{Any, Any})
    mesh_struct = params["mesh_struct"]
    x, y, z = mesh_struct.get_coords()
    num_nodes = length(x)
    reference = Matrix{MTScalar}(undef, 3, num_nodes)
    temperature = Vector{MTScalar}(undef, num_nodes)
    rate = Vector{MTScalar}(undef, num_nodes)
    for node ∈ 1 : num_nodes
        reference[:, node] = [x[node], y[node], z[node]]
        temperature[node] = 0.0
        rate[node] = 0.0
    end
    materials_file = params["material"]
    material_params = YAML.load_file(materials_file)
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    if (num_blks_params ≠ num_blks)
        error("number of blocks in mesh ", params["mesh"], " (", num_blks,
        ") must be equal to number of blocks in materials file ", params["material"],
        " (", num_blks_params, ")")
    end
    elem_blk_names = mesh_struct.get_elem_blk_names()
    materials = Vector{Thermal}(undef, 0)
    for elem_blk_name ∈ elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        push!(materials, material_model)
    end
    time = 0.0
    nodal_dofs = [free::DOF for _ ∈ 1 : num_nodes]
    DynamicHeat(params, materials, reference, temperature, rate, nodal_dofs, time)
end

function create_model(params::Dict{Any, Any})
    model_name = params["model"]
    if model_name == "static solid"
        return StaticSolid(params)
    elseif model_name == "dynamic solid"
        return DynamicSolid(params)
    elseif model_name == "static heat"
        return StaticHeat(params)
    elseif model_name == "dynamic heat"
        return DynamicHeat(params)
    else
        error("Unknown type of model : ", model_name)
    end
end

function potential_energy(model::SolidMechanics)
    params = model.params
    materials = model.materials
    mesh_struct = params["mesh_struct"]
    body_energy = 0.0
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1 : num_blks
        material = materials[blk_index]
        blk_id = elem_blk_ids[blk_index]
        elem_type = mesh_struct.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        _, dNdξ, elem_weights = isoparametric(elem_type, num_points)
        blk_conn = mesh_struct.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        num_elem_nodes = blk_conn[3]
        for blk_elem_index ∈ 1 : num_blk_elems
            node_indices = (blk_elem_index - 1) * num_elem_nodes + 1 : blk_elem_index * num_elem_nodes 
            elem_ref_pos = model.reference[:, node_indices]
            elem_cur_pos = model.current[:, node_indices]
            element_energy = 0.0
            for point ∈ 1 : num_points
                dXdξ = dNdξ[:, :, point] * elem_ref_pos'
                dxdξ = dNdξ[:, :, point] * elem_cur_pos'
                dxdX = dXdξ \ dxdξ
                j = det(dXdξ)
                F = MTTensor(dxdX)
                W, _, _ = constitutive(material, F)
                w = elem_weights[point]
                element_energy += W * j * w
            end
            body_energy += element_energy
        end
    end
    return body_energy
end

function energy_force_stiffness(model::SolidMechanics)
    params = model.params
    materials = model.materials
    mesh_struct = params["mesh_struct"]
    x, _, _ = mesh_struct.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    total_energy = 0.0
    total_internal_force = spzeros(num_dof)
    total_stiffness = spzeros(num_dof, num_dof)
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1 : num_blks
        material = materials[blk_index]
        blk_id = elem_blk_ids[blk_index]
        elem_type = mesh_struct.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        _, dNdξ, elem_weights = isoparametric(elem_type, num_points)
        blk_conn = mesh_struct.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        num_elem_nodes = blk_conn[3]
        for blk_elem_index ∈ 1 : num_blk_elems
            node_indices = (blk_elem_index - 1) * num_elem_nodes + 1 : blk_elem_index * num_elem_nodes 
            elem_ref_pos = model.reference[:, node_indices]
            elem_cur_pos = model.current[:, node_indices]
            element_energy = 0.0
            element_internal_force = spzeros(3 * num_elem_nodes)
            element_stiffness = spzeros(3 * num_elem_nodes, 3 * num_elem_nodes)
            elem_dofs = zeros(Int64, 3 * num_elem_nodes)
            elem_dofs[1 : num_elem_nodes] = 3 .* node_indices .- 2
            elem_dofs[num_elem_nodes + 1 : 2 * num_elem_nodes] = 3 .* node_indices .- 1
            elem_dofs[2 * num_elem_nodes + 1 : 3 * num_elem_nodes] = 3 .* node_indices
            for point ∈ 1 : num_points
                dXdξ = dNdξ[:, :, point] * elem_ref_pos'
                dxdξ = dNdξ[:, :, point] * elem_cur_pos'
                dxdX = dXdξ \ dxdξ
                dNdX = dXdξ \ dNdξ[:, :, point]
                B = gradient_operator(dNdX)
                j = det(dXdξ)
                F = MTTensor(dxdX)
                W, P, A = constitutive(material, F)
                stress = reshape(P', 9, 1)
                moduli = second_from_fourth(A)
                w = elem_weights[point]
                element_energy += W * j * w
                element_internal_force += B' * stress * j * w
                element_stiffness += B' * moduli * B * j * w
            end
            total_energy += element_energy
            total_internal_force[elem_dofs] += element_internal_force
            total_stiffness[elem_dofs, elem_dofs] += element_stiffness
        end
    end
    return total_energy, total_internal_force, total_stiffness
end

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
    model.nodal_dofs = nodal_dofs
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
    model.nodal_dofs = nodal_dofs
end