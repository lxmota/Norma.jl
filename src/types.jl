struct ProblemParams
    problem_type::String
    num_domains::Int
    num_steps::Int
    start_time::Float64
    stop_time::Float64
    rel_tol_schwarz::Float64
    abs_tol_schwarz::Float64
    rel_tol_newton::Float64
    abs_tol_newton::Float64
    max_iter_schwarz::Int
    max_iter_newton::Int
    step_interval::Int
    schwarz_interval::Int
    newton_interval::Int
end

struct BoundaryCondition
    nodesets::Vector{Int}
    dof::Vector{Int}
    time::Vector{Float64}
    value::Vector{Float64}
end

struct Mesh
    name::String
    dimension::Int
    num_nodes::Int
    num_elements::Int
    element_type::String
    nodes_per_element::Int
    num_int::Int
    nodesets::Vector{Vector{Int}}
    dbc::BoundaryCondition
    nbc::BoundaryCondition
    sbc::BoundaryCondition
    coordinates::Array{Float64}
    connectivity::Array{Int}
    function Mesh(name, dimension, num_nodes, num_elements, element_type,
                  nodes_per_element, num_int)
        nodesets = Vector{Int}[]
        dbc = BoundaryCondition([], [], [], [])
        nbc = BoundaryCondition([], [], [], [])
        sbc = BoundaryCondition([], [], [], [])
        coordinates = zeros(dimension, num_nodes)
        connectivity = zeros(num_elements, nodes_per_element)
        new(name, dimension, num_nodes, num_elements, element_type,
            nodes_per_element, num_int, nodesets, dbc, nbc, sbc,
            coordinates, connectivity)
    end
end

struct Domain
    Na::Array{Float64}
    dNadξ::Array{Float64}
    w::Array{Float64}
    function Domain(mesh::Mesh)
        Na, dNadξ, w = isoparametric(mesh.element_type, mesh.num_int)        
    end
end
