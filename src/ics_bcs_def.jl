abstract type BoundaryCondition end
abstract type SchwarzBoundaryCondition <: BoundaryCondition end
abstract type RegularBoundaryCondition <: BoundaryCondition end
abstract type ContactSchwarzBoundaryCondition <: SchwarzBoundaryCondition end
abstract type RegularSchwarzBoundaryCondition <: SchwarzBoundaryCondition end
abstract type InitialCondition end

using Symbolics

mutable struct SMDirichletBC <: RegularBoundaryCondition
    node_set_name::String
    offset::Int64
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_num::Num
    velo_num::Num
    acce_num::Num
end

mutable struct SMNeumannBC <: RegularBoundaryCondition
    side_set_name::String
    offset::Int64
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    traction_num::Num
end

mutable struct SMContactSchwarzBC <: ContactSchwarzBoundaryCondition
    side_set_name::String
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    coupled_subsim::Simulation
    coupled_bc_index::Int64
    coupled_mesh::PyObject
    coupled_block_id::Int64
    coupled_side_set_id::Int64
    is_dirichlet::Bool
    projection_operator::Matrix{Float64}
end

mutable struct SMSchwarzDBC <: RegularSchwarzBoundaryCondition
    node_set_name::String
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    coupled_subsim::Simulation
    coupled_mesh::PyObject
    coupled_block_id::Int64
    coupled_nodes_indices::Vector{Vector{Int64}}
    interpolation_function_values::Vector{Vector{Float64}}
end