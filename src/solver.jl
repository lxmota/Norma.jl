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
    gradient::SparseVector{Float64, Int64}
    hessian::SparseMatrixCSC{Float64, Int64}
    initial_guess::SparseVector{Float64, Int64}
    initial_norm::Float64
    converged::Bool
    failed::Bool
    step::Step
end

struct NewtonStep <: Step
end

function create_step(solver_params::Dict{Any, Any})
    step_name = solver_params["step"]
    if step_name == "full Newton"
        return NewtonStep()
    else
        error("Unknown type of solver step: ", step_name)
    end
end

function HessianMinimizer(params::Dict{Any, Any})
    solver_params = params["solver"]
    mesh_struct = params["mesh_struct"]
    x, _, _ = mesh_struct.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    minimum_iterations = solver_params["minimum iterations"]
    maximum_iterations = solver_params["maximum iterations"]
    iteration_number = 1
    absolute_tolerance = solver_params["absolute tolerance"]
    relative_tolerance = solver_params["relative tolerance"]
    absolute_error = 0.0
    relative_error = 0.0
    value = 0.0
    gradient = spzeros(num_dof)
    hessian = spzeros(num_dof, num_dof)
    initial_guess = spzeros(num_dof)
    initial_norm = 0.0
    converged = false
    failed = false
    step = create_step(solver_params)
    HessianMinimizer(minimum_iterations, maximum_iterations, iteration_number, absolute_tolerance,
    relative_tolerance, absolute_error, relative_error, value, gradient, hessian, initial_guess,
    initial_norm, converged, failed, step)
   end

function create_solver(params::Dict{Any, Any})
    solver_params = params["solver"]
    solver_name = solver_params["type"]
    if solver_name == "Hessian minimizer"
        return HessianMinimizer(params)
    else
        error("Unknown type of solver : ", solver_name)
    end
end

function evaluate(solver::HessianMinimizer, model::SolidMechanics, solution::SparseVector{Float64, Int64})
    num_nodes = length(solution) ÷ 3
    for node ∈ 1 : num_nodes
        model.current[:, node] = [solution[3 * node - 2], solution[3 * node - 1], solution[3 * node]]
    end
    solver.value, solver.gradient, solver.hessian = energy_force_stiffness(model)
end

function compute_step(solver::HessianMinimizer, model::SolidMechanics, step_type::NewtonStep, solution::SparseVector{Float64, Int64})
    evaluate(solver, model, solution)
    step = - solver.hessian \ solver.gradient
    return step
end

function update_convergence_criterion(solver::HessianMinimizer, absolute_error::Float64)
    solver.absolute_error = absolute_error
    solver.relative_error = solver.initial_norm > 0.0 ? absolute_error / solver.initial_norm : 0.0
    converged_absolute = solver.absolute_error ≤ solver.absolute_tolerance
    converged_relative = solver.relative_error ≤ solver.relative_tolerance
    solver.converged = converged_absolute || converged_relative
end

function continue_solve(solver::HessianMinimizer)
    if solver.failed == true
        return false
    end
    zero_residual = solver.absolute_error == 0.0
    if zero_residual == true
        return false
    end
    exceeds_minimum_iterations = solver.iteration_number > solver.minimum_iterations
    if exceeds_minimum_iterations == false
        return true
    end
    exceeds_maximum_iterations = solver.iteration_number > solver.maximum_iterations
    if exceeds_maximum_iterations == true
        return false
    end
    continue_solving = solver.converged == false
    return continue_solving
end

function solve(model::SolidMechanics, solver::HessianMinimizer, solution::SparseVector{Float64, Int64})
    solver.initial_guess = solution
    evaluate(solver, model, solution)
    residual = solver.gradient
    solver.initial_norm = norm(residual)
    solver.iteration_number = 1
    solver.failed = solver.failed || model.failed
    step_type = solver.step
    update_convergence_criterion(solver, solver.initial_norm)
    while continue_solve(solver) == true
        step = compute_step(solver, model, step_type, solution)
        solution += step
        evaluate(solver, model, solution)
        residual = solver.gradient
        norm_residual = norm(residual)
        update_convergence_criterion(solver, norm_residual)
        solver.iteration_number += 1
    end
end