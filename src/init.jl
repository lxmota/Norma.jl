include("model.jl")
include("solver.jl")
include("time.jl")

function init(params::Dict{Any, Any})
    model = create_model(params)
    params["model_struct"] = model
    solver = create_solver(params)
    params["solver_struct"] = solver
    time_integrator = create_time_integrator(params)
    params["time_integrator_struct"] = time_integrator
end