abstract type Model end
using SparseArrays

@enum Kinematics begin
    Undefined
    Infinitesimal
    Finite
end

mutable struct SolidMechanics <: Model
    mesh::ExodusDatabase
    materials::Vector{Solid}
    reference::Matrix{Float64}
    current::Matrix{Float64}
    velocity::Matrix{Float64}
    acceleration::Matrix{Float64}
    internal_force::Vector{Float64}
    boundary_force::Vector{Float64}
    boundary_conditions::Vector{BoundaryCondition}
    stress::Vector{Vector{Vector{Vector{Float64}}}}
    stored_energy::Vector{Vector{Float64}}
    free_dofs::BitVector
    time::Float64
    failed::Bool
    mesh_smoothing::Bool
    smooth_reference::String
    inclined_support::Bool
    global_transform::SparseArrays.Matrix{Float64}
    kinematics::Kinematics
end

# TODO: Add potential energy as in the above
mutable struct HeatConduction <: Model
    mesh::ExodusDatabase
    materials::Vector{Vector}
    reference::Matrix{Float64}
    temperature::Vector{Float64}
    rate::Vector{Float64}
    internal_heat_flux::Vector{Float64}
    boundary_heat_flux::Vector{Float64}
    boundary_conditions::Vector{BoundaryCondition}
    flux::Vector{Vector{Vector{Vector{Float64}}}}
    stored_energy::Vector{Vector{Float64}}
    free_dofs::BitVector
    time::Float64
    failed::Bool
end
