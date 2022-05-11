include("constitutive.jl")
include("field.jl")

abstract type Model end
abstract type SolidMechanics <: Model end
abstract type HeatConduction <: Model end

mutable struct StaticSolid <: SolidMechanics
    params::Dict{Any, Any}
    materials::Vector{Material}
    reference::VectorField
    current::VectorField
end

function StaticSolid(params::Dict{Any, Any})
    mesh_struct = params["mesh_struct"]
    x, y, z = mesh_struct.get_coords()
    num_nodes = length(x)
    reference_vector = Vector{MTVector}(undef, num_nodes)
    current_vector = Vector{MTVector}(undef, num_nodes)
    for node ∈ 1 : num_nodes
        reference_vector[node] = MTVector(x[node], y[node], z[node])
        current_vector[node] = MTVector(x[node], y[node], z[node])
    end
    reference = VectorField("reference configuration", reference_vector)
    current = VectorField("current configuration", current_vector)
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
    materials = Vector{Material}(undef, 0)
    for elem_blk_name ∈ elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        push!(materials, material_model)
    end
    StaticSolid(params, materials, reference, current)
end

mutable struct DynamicSolid <: SolidMechanics
    params::Dict{Any, Any}
    materials::Vector{Material}
    reference::VectorField
    current::VectorField
    velocity::VectorField
    acceleration::VectorField
end

mutable struct StaticHeat <: HeatConduction
    params::Dict{Any, Any}
    reference::VectorField
    temperature::ScalarField
end

mutable struct DynamicHeat <: HeatConduction
    params::Dict{Any, Any}
    reference::VectorField
    temperature::ScalarField
    rate::ScalarField
end

function create_model(params::Dict{Any, Any})
    model_name = params["model"]
    mesh_struct = params["mesh_struct"]
    x, y, z = mesh_struct.get_coords()
    num_nodes = length(x)
    reference = Vector{MTVector}(undef, num_nodes)
    for node ∈ 1 : num_nodes
        reference[node] = MTVector(x[node], y[node], z[node])
    end
    reference_field = VectorField("reference configuration", reference)

    if model_name == "static solid" || model_name == "dynamic solid"
        current = Vector{MTVector}(undef, num_nodes)
        for node ∈ 1 : num_nodes
            current[node] = MTVector(x[node], y[node], z[node])
        end
        current_field = VectorField("current configuration", current)
        materials_file = params["material"]
        material_params = YAML.load_file(materials_file)
    end
    if model_name == "static heat" || model_name == "dynamic heat"
        temperature = Vector{MTScalar}(undef, num_nodes)
        for node ∈ 1 : num_nodes
            temperature[node] = 0.0
        end
        temperature_field = ScalarField("temperature", temperature)
    end
  
    if model_name == "static solid"
        return StaticSolid(params)
    elseif model_name == "dynamic solid"
        velocity = Vector{MTVector}(undef, num_nodes)
        acceleration = Vector{MTVector}(undef, num_nodes)
        for node ∈ 1 : num_nodes
            velocity[node] = acceleration[node] = MTVector(0.0, 0.0, 0.0)
        end
        velocity_field = VectorField("velocity", velocity)
        acceleration_field = VectorField("acceleration", acceleration)
        return DynamicSolid(params, material_params, reference_field, current_field,
            velocity_field, acceleration_field)
    elseif model_name == "static heat"
        return StaticHeat(params, reference_field, temperature_field)
    elseif model_name == "dynamic heat"
        rate = Vector{MTScalar}(undef, num_nodes)
        for node ∈ 1 : num_nodes
            rate[node] = 0.0
        end
        rate_field = ScalarField("temperature rate", rate)
        return DynamicHeat(params, reference_field, temperature_field, rate_field)
    else
        error("Unknown type of model : ", model_name)
    end
end

function potential_energy(model::SolidMechanics)
    params = model.params
    material_params = model.material_params
    mesh_struct = params["mesh_struct"]
    total = 0.0
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1 : num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = mesh_struct.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        N, dNdξ, w = isoparametric(elem_type, num_points)
        blk_conn = mesh_struct.get_elem_connectivity(blk_id)
        elem_conn = blk_conn[1]
        num_blk_elems = blk_conn[2]
        num_elem_nodes = blk_conn[3]
        for blk_elem_index ∈ 1 : num_blk_elems
            elem_nodes = (blk_elem_index - 1) * num_elem_nodes + 1 : blk_elem_index * num_elem_nodes 
            connectvity = elem_conn[elem_nodes]
            elem_ref_pos = model.reference.value[elem_nodes]
            elem_cur_pos = model.current.value[elem_nodes]
            for point ∈ 1 : num_points
            end
        end
    end
end