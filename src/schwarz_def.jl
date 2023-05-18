abstract type SchwarzController end

mutable struct SolidSchwarzController <: SchwarzController
    num_domains::Int64
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
    prev_time::Float64
    stop::Int64
    converged::Bool
    stop_disp::Vector{Vector{Float64}}
    stop_velo::Vector{Vector{Float64}}
    stop_acce::Vector{Vector{Float64}}
    stop_traction_force::Vector{Vector{Float64}}
    schwarz_disp::Vector{Vector{Float64}}
    schwarz_velo::Vector{Vector{Float64}}
    schwarz_acce::Vector{Vector{Float64}}
    time_hist::Vector{Vector{Float64}}
    disp_hist::Vector{Vector{Vector{Float64}}}
    velo_hist::Vector{Vector{Vector{Float64}}}
    acce_hist::Vector{Vector{Vector{Float64}}}
    traction_force_hist::Vector{Vector{Vector{Float64}}}
    schwarz_contact::Bool
    active_contact::Bool
end