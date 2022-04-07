include("model.jl")

function init(params::Dict{Any, Any})
end

function init(model::StaticSolid, params::Dict{Any, Any})
    mesh_struct = params["mesh_struct"]
end