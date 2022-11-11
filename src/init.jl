include("model.jl")
include("solver.jl")

function init(params::Dict{Any, Any})
    model = create_model(params)
    params["model_struct"] = model
    solver = create_solver(params)
    params["solver_struct"] = solver
end