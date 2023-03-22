include("constitutive_def.jl")
include("exodus_def.jl")
include("ics_bcs_def.jl")
include("model_def.jl")
include("time_def.jl")
include("solver_def.jl")
include("schwarz_def.jl")

abstract type Simulation end

struct SingleDomainSimulation <: Simulation
    name::String
    params::Dict{Any,Any}
    integrator::TimeIntegrator
    solver::Solver
    model::Model
end

struct MultiDomainSimulation <: Simulation
    name::String
    params::Dict{Any,Any}
    schwarz_controller::SchwarzController
    subsims::Vector{SingleDomainSimulation}
    subsim_name_index_map::Dict{String,Int64}
end