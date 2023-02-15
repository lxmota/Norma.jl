abstract type Solver end
abstract type Minimizer <: Solver end
abstract type Step end

mutable struct HessianMinimizer <: Minimizer
    minimum_iterations::Int64
    maximum_iterations::Int64
    iteration_number::Int64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    absolute_error::Float64
    relative_error::Float64
    value::Float64
    gradient::Vector{Float64}
    hessian::SparseMatrixCSC{Float64, Int64}
    solution::Vector{Float64}
    free_dofs::BitVector
    initial_norm::Float64
    converged::Bool
    failed::Bool
    step::Step
end

struct NewtonStep <: Step
end