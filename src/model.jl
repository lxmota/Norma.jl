include("field.jl")

abstract type Model end
abstract type SolidMechanics <: Model end
abstract type HeatConduction <: Model end

mutable struct StaticSolid <: SolidMechanics
    reference::VectorField
    current::VectorField
end

mutable struct DynamicSolid <: SolidMechanics
    reference::VectorField
    current::VectorField
    velocity::VectorField
    acceleration::VectorField
end

mutable struct StaticHeat <: HeatConduction
    reference::VectorField
    temperature::ScalarField
end

mutable struct DynamicHeat <: HeatConduction
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
    end
    if model_name == "static heat" || model_name == "dynamic heat"
        temperature = Vector{MTScalar}(undef, num_nodes)
        for node ∈ 1 : num_nodes
            temperature[node] = 0.0
        end
        temperature_field = ScalarField("temperature", temperature)
    end
  
    if model_name == "static solid"
        return StaticSolid(reference_field, current_field)
    elseif model_name == "dynamic solid"
        velocity = Vector{MTVector}(undef, num_nodes)
        acceleration = Vector{MTVector}(undef, num_nodes)
        for node ∈ 1 : num_nodes
            velocity[node] = acceleration[node] = MTVector(0.0, 0.0, 0.0)
        end
        velocity_field = VectorField("velocity", velocity)
        acceleration_field = VectorField("acceleration", acceleration)
        return DynamicSolid(reference_field, current_field, velocity_field, acceleration_field)
    elseif model_name == "static heat"
        return StaticHeat(reference_field, temperature_field)
    elseif model_name == "dynamic heat"
        rate = Vector{MTScalar}(undef, num_nodes)
        for node ∈ 1 : num_nodes
            rate[node] = 0.0
        end
        rate_field = ScalarField("temperature rate", rate)
        return DynamicHeat(reference_field, temperature_field, rate_field)
    else
        error("Unknown type of model ", model_name)
    end
end