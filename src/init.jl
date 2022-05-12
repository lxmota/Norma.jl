include("model.jl")

function init(params::Dict{Any, Any})
    model = create_model(params)
    params["model_struct"] = model
end