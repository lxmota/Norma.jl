abstract type TimeIntegrator end

mutable struct QuasiStatic <: TimeIntegrator
    initial_time::Float64
    final_time::Float64
    time_step::Float64
end

mutable struct Newmark <: TimeIntegrator
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    β::Float64
    γ::Float64
end

function create_time_integrator(params::Dict{Any, Any})
    solver_params = params["time integrator"]
    solver_name = solver_params["type"]
    if solver_name == "Hessian minimizer"
        return HessianMinimizer(params)
    else
        error("Unknown type of solver : ", solver_name)
    end
end
