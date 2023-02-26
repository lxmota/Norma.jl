abstract type SchwarzController end

mutable struct StaticSchwarzController <: SchwarzController
    domains::Vector{String}
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    previous_solution::Vector{Float64}
end

mutable struct DynamicSchwarzController <: SchwarzController
    domains::Vector{String}
    initial_time::Float64
    final_time::Float64
    time_step::Float64
    previous_solution::Vector{Float64}
    previous_sol_dot::Vector{Float64}
    previous_sol_dotdot::Vector{Float64}
end