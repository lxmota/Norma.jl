abstract type SchwarzController end

mutable struct StaticSchwarzController <: SchwarzController
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    time::Float64
    stop::Int64
    previous_solutions::Vector{Vector{Float64}}
end

mutable struct DynamicSchwarzController <: SchwarzController
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    time::Float64
    stop::Int64
    previous_solutions::Vector{Vector{Float64}}
    previous_sol_dots::Vector{Vector{Float64}}
    previous_sol_dotdots::Vector{Vector{Float64}}
end