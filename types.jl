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

struct Mesh
    name
    dimension
    num_nodes
    num_elements
    element_type
    nodes_per_element
    num_int
    coordinates
    connectivity
    nodesets
    dbc
    nbc
    sbc
end

struct Domain
    Na
    dNadξ
    w
end
