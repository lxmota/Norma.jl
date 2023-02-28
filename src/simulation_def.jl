include("constitutive_def.jl")
include("exodus_def.jl")
include("model_def.jl")
include("time_def.jl")
include("solver_def.jl")
include("schwarz_def.jl")

abstract type Simulation end

struct SingleDomainSimulation <: Simulation
    params::Dict{Any,Any}
    integrator::TimeIntegrator
    solver::Solver
    model::Model
end

struct MultiDomainSimulation <: Simulation
    params::Dict{Any,Any}
    schwarz_controller::SchwarzController
    sub_simulations::Vector{SingleDomainSimulation}
end