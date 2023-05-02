abstract type BoundaryCondition end
abstract type InitialCondition end

using Symbolics

mutable struct SMDirichletBC <: BoundaryCondition
    node_set_name::String
    offset::Int64
    node_set_id::Int64
    node_set_node_indices::Vector{Int64}
    disp_num::Num
    velo_num::Num
    acce_num::Num
end

mutable struct SMNeumannBC <: BoundaryCondition
    side_set_name::String
    offset::Int64
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    traction_num::Num
end

mutable struct SMSchwarzContactBC <: BoundaryCondition
    side_set_name::String
    side_set_id::Int64
    num_nodes_per_side::Vector{Int64}
    side_set_node_indices::Vector{Int64}
    coupled_subsim::Simulation
    coupled_mesh::PyObject
    coupled_block_id::Int64
    coupled_side_set_id::Int64
    is_dirichlet::Bool
end