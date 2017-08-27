struct ProblemParams
    problem_type
    num_domains
    num_steps
    start_time
    stop_time
    rel_tol_schwarz
    abs_tol_schwarz
    rel_tol_newton
    abs_tol_newton
    max_iter_schwarz
    max_iter_newton
    step_interval
    schwarz_interval
    newton_interval
    dimension
    dof
end

struct Domain
    name
    dimension
    num_nodes
    num_elements
    nodes_per_element
    num_int
    element_type
    mesh
    dbc
    nbc
    sbc
end

struct Mesh
    coordinates
    connectivity
    element_type
    nodesets
end
