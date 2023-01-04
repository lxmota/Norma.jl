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
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
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
    gradient = zeros(num_dof)
    hessian = spzeros(num_dof, num_dof)
    free_dofs = trues(num_dof)
    initial_guess = zeros(num_dof)
    initial_norm = 0.0
    converged = false
    failed = false
    step = create_step(solver_params)
    HessianMinimizer(minimum_iterations, maximum_iterations, iteration_number,
    absolute_tolerance, relative_tolerance, absolute_error, relative_error,
    value, gradient, hessian, initial_guess, free_dofs,
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

function copy_solution_source_targets(solver::HessianMinimizer, integrator::QuasiStatic, model::SolidMechanics)
    displacement = solver.solution
    integrator.displacement = displacement
    _, num_nodes = size(model.reference)
    for node ∈ 1 : num_nodes
        nodal_displacement = displacement[3 * node - 2 : 3 * node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
end

function copy_solution_source_targets(model::SolidMechanics, integrator::QuasiStatic, solver::HessianMinimizer)
    _, num_nodes = size(model.reference)
    for node ∈ 1 : num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        integrator.displacement[3 * node - 2 : 3 * node] = nodal_displacement
    end
    solver.solution = integrator.displacement
end

function copy_solution_source_targets(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    displacement = integrator.displacement
    solver.solution = displacement
    _, num_nodes = size(model.reference)
    for node ∈ 1 : num_nodes
        nodal_displacement = displacement[3 * node - 2 : 3 * node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
end

function copy_solution_source_targets(solver::HessianMinimizer, integrator::Newmark, model::SolidMechanics)
    displacement = solver.solution
    integrator.displacement = displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    _, num_nodes = size(model.reference)
    for node ∈ 1 : num_nodes
        nodal_displacement = displacement[3 * node - 2 : 3 * node]
        nodal_velocity = velocity[3 * node - 2 : 3 * node]
        nodal_acceleration = acceleration[3 * node - 2 : 3 * node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(model::SolidMechanics, integrator::Newmark, solver::HessianMinimizer)
    _, num_nodes = size(model.reference)
    for node ∈ 1 : num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        nodal_velocity = model.velocity[:, node]
        nodal_acceleration = model.acceleration[:, node]
        integrator.displacement[3 * node - 2 : 3 * node] = nodal_displacement
        integrator.velocity[3 * node - 2 : 3 * node] = nodal_velocity
        integrator.acceleration[3 * node - 2 : 3 * node] = nodal_acceleration
    end
    solver.solution = integrator.displacement
end

function copy_solution_source_targets(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    solver.solution = displacement
    _, num_nodes = size(model.reference)
    for node ∈ 1 : num_nodes
        nodal_displacement = displacement[3 * node - 2 : 3 * node]
        nodal_velocity = velocity[3 * node - 2 : 3 * node]
        nodal_acceleration = acceleration[3 * node - 2 : 3 * node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function evaluate(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    strain_energy, internal_force, external_force, stiffness_matrix, _ = evaluate(model)
    solver.value = strain_energy
    solver.gradient = internal_force - external_force
    solver.hessian = stiffness_matrix
end

function evaluate(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    strain_energy, internal_force, external_force, stiffness_matrix, mass_matrix = evaluate(model)
    β = integrator.β
    Δt = integrator.time_step
    implicit = β > 0.0
    if implicit == true
        inertial_force = mass_matrix * integrator.acceleration
        kinetic_energy = 0.5 * integrator.velocity' * mass_matrix * integrator.velocity
        solver.hessian = stiffness_matrix + mass_matrix / β / Δt / Δt
    else
        mass_vector = sum(mass_matrix, dim = 2)
        inertial_force = mass_vector .* integrator.acceleration
        kinetic_energy = 0.5 * mass_vector .* integrator.velocity .* integrator.velocity
    end
    solver.value = strain_energy - external_force' * integrator.displacement + kinetic_energy
    solver.gradient = internal_force - external_force + inertial_force
end

function compute_step(solver::HessianMinimizer, step_type::NewtonStep)
    step = zeros(length(solver.gradient))
    free_dofs = solver.free_dofs
    step[free_dofs] = - solver.hessian[free_dofs, free_dofs] \ solver.gradient[free_dofs]
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

function solve(integrator::Any, model::SolidMechanics, solver::HessianMinimizer)
    copy_solution_source_targets(model, integrator, solver)
    predict(integrator, solver, model)
    evaluate(integrator, solver, model)
    residual = solver.gradient
    norm_residual = norm(residual[solver.free_dofs])
    solver.initial_norm = norm_residual
    solver.iteration_number = 1
    solver.failed = solver.failed || model.failed
    step_type = solver.step
    update_convergence_criterion(solver, solver.initial_norm)
    println("initial |R|=", norm_residual, ", |Δ|=", norm(solver.solution))
    while true
        step = compute_step(solver, step_type)
        solver.solution += step
        copy_solution_source_targets(solver, integrator, model)
        correct(integrator, solver, model)
        evaluate(integrator, solver, model)
        residual = solver.gradient
        norm_step = norm(step)
        norm_residual = norm(residual[solver.free_dofs])
        println("iter=", solver.iteration_number, ", |R|=", norm_residual, ", |Δ|=", norm(solver.solution), ", |δΔ|=", norm_step)
        update_convergence_criterion(solver, norm_residual)
        solver.iteration_number += 1
        if continue_solve(solver) == false
            break
        end
    end
end