include("constitutive.jl")

abstract type Model end
abstract type SolidMechanics <: Model end
abstract type HeatConduction <: Model end

mutable struct StaticSolid <: SolidMechanics
    params::Dict{Any, Any}
    materials::Vector{Solid}
    reference::Matrix{MTScalar}
    current::Matrix{MTScalar}
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
    StaticSolid(params, materials, reference, current, time)
end

mutable struct DynamicSolid <: SolidMechanics
    params::Dict{Any, Any}
    materials::Vector{Solid}
    reference::Matrix{MTScalar}
    current::Matrix{MTScalar}
    velocity::Matrix{MTScalar}
    acceleration::Matrix{MTScalar}
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
    DynamicSolid(params, materials, reference, current, velocity, acceleration, time)
end

mutable struct StaticHeat <: HeatConduction
    params::Dict{Any, Any}
    materials::Vector{Thermal}
    reference::Matrix{MTScalar}
    temperature::Vector{MTScalar}
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
    StaticHeat(params, materials, reference, temperature, time)
end

mutable struct DynamicHeat <: HeatConduction
    params::Dict{Any, Any}
    materials::Vector{Vector}
    reference::Matrix{MTScalar}
    temperature::Vector{MTScalar}
    rate::Vector{MTScalar}
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
    DynamicHeat(params, materials, reference, temperature, rate, time)
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
        elem_conn = blk_conn[1]
        num_blk_elems = blk_conn[2]
        num_elem_nodes = blk_conn[3]
        for blk_elem_index ∈ 1 : num_blk_elems
            node_indices = (blk_elem_index - 1) * num_elem_nodes + 1 : blk_elem_index * num_elem_nodes 
            elem_nodes = elem_conn[node_indices]
            elem_ref_pos = model.reference[:, elem_nodes]
            elem_cur_pos = model.current[:, elem_nodes]
            element_energy = 0.0
            for point ∈ 1 : num_points
                dXdξ = dNdξ[:, :, point] * transpose(elem_ref_pos)
                dxdξ = dNdξ[:, :, point] * transpose(elem_cur_pos)
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

function potential_energy(position, model::SolidMechanics)
    mesh_struct = model.params["mesh_struct"]
    num_nodes = mesh_struct.num_nodes()
    for node ∈ 1 : num_nodes
        model.current[1, node] = position[3 * node - 2]
        model.current[2, node] = position[3 * node - 1]
        model.current[3, node] = position[3 * node]
    end
    return potential_energy(model)
end