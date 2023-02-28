abstract type Model end

mutable struct SolidMechanics <: Model
    mesh::PyObject
    materials::Vector{Solid}
    reference::Matrix{Float64}
    current::Matrix{Float64}
    velocity::Matrix{Float64}
    acceleration::Matrix{Float64}
    boundary_tractions_force::Vector{Float64}
    stress::Vector{Vector{Vector{Vector{Float64}}}}
    free_dofs::BitVector
    time::Float64
    failed::Bool
end

mutable struct HeatConduction <: Model
    mesh::PyObject
    materials::Vector{Vector}
    reference::Matrix{Float64}
    temperature::Vector{Float64}
    rate::Vector{Float64}
    boundary_heat_flux::Vector{Float64}
    flux::Vector{Vector{Vector{Vector{Float64}}}}
    free_dofs::BitVector
    time::Float64
    failed::Bool
end