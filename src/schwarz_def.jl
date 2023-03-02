abstract type SchwarzController end

mutable struct SolidStaticSchwarzController <: SchwarzController
    minimum_iterations::Int64
    maximum_iterations::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    absolute_error::Float64
    relative_error::Float64
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    time::Float64
    stop::Int64
    converged::Bool
    prev_disp::Vector{Vector{Float64}}
end

mutable struct SolidDynamicSchwarzController <: SchwarzController
    minimum_iterations::Int64
    maximum_iterations::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    absolute_error::Float64
    relative_error::Float64
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    time::Float64
    stop::Int64
    converged::Bool
    prev_disp::Vector{Vector{Float64}}
    prev_velo::Vector{Vector{Float64}}
    prev_acce::Vector{Vector{Float64}}
end