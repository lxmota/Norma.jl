function create_step(solver_params::Dict{Any,Any})
    step_name = solver_params["step"]
    if step_name == "full Newton"
        return NewtonStep()
    elseif step_name == "explicit"
        return ExplicitStep()
    else
        error("Unknown type of solver step: ", step_name)
    end
end

function HessianMinimizer(params::Dict{Any,Any})
    solver_params = params["solver"]
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    minimum_iterations = solver_params["minimum iterations"]
    maximum_iterations = solver_params["maximum iterations"]
    absolute_tolerance = solver_params["absolute tolerance"]
    relative_tolerance = solver_params["relative tolerance"]
    absolute_error = 0.0
    relative_error = 0.0
    value = 0.0
    gradient = zeros(num_dof)
    hessian = spzeros(num_dof, num_dof)
    initial_guess = zeros(num_dof)
    initial_norm = 0.0
    converged = false
    failed = false
    step = create_step(solver_params)
    HessianMinimizer(minimum_iterations, maximum_iterations,
        absolute_tolerance, relative_tolerance, absolute_error, relative_error,
        value, gradient, hessian, initial_guess,
        initial_norm, converged, failed, step)
end

function ExplicitSolver(params::Dict{Any,Any})
    solver_params = params["solver"]
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    value = 0.0
    gradient = zeros(num_dof)
    lumped_hessian = zeros(num_dof)
    initial_guess = zeros(num_dof)
    initial_norm = 0.0
    converged = false
    failed = false
    step = create_step(solver_params)
    ExplicitSolver(value, gradient, lumped_hessian, initial_guess,
        initial_norm, converged, failed, step)
end

function create_solver(params::Dict{Any,Any})
    solver_params = params["solver"]
    solver_name = solver_params["type"]
    if solver_name == "Hessian minimizer"
        return HessianMinimizer(params)
    elseif solver_name == "explicit solver"
        return ExplicitSolver(params)
    else
        error("Unknown type of solver : ", solver_name)
    end
end

function copy_solution_source_targets(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    displacement = integrator.displacement
    solver.solution = displacement
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
end

function copy_solution_source_targets(solver::HessianMinimizer, model::SolidMechanics, integrator::QuasiStatic)
    displacement = solver.solution
    integrator.displacement = displacement
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
    end
end

function copy_solution_source_targets(model::SolidMechanics, integrator::QuasiStatic, solver::HessianMinimizer)
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        integrator.displacement[3*node-2:3*node] = nodal_displacement
    end
    solver.solution = integrator.displacement
end

function copy_solution_source_targets(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    solver.solution = displacement
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        nodal_velocity = velocity[3*node-2:3*node]
        nodal_acceleration = acceleration[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(solver::HessianMinimizer, model::SolidMechanics, integrator::Newmark)
    displacement = solver.solution
    integrator.displacement = displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        nodal_velocity = velocity[3*node-2:3*node]
        nodal_acceleration = acceleration[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(model::SolidMechanics, integrator::Newmark, solver::HessianMinimizer)
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        nodal_velocity = model.velocity[:, node]
        nodal_acceleration = model.acceleration[:, node]
        integrator.displacement[3*node-2:3*node] = nodal_displacement
        integrator.velocity[3*node-2:3*node] = nodal_velocity
        integrator.acceleration[3*node-2:3*node] = nodal_acceleration
    end
    solver.solution = integrator.displacement
end

function copy_solution_source_targets(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = integrator.acceleration
    solver.solution = acceleration
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        nodal_velocity = velocity[3*node-2:3*node]
        nodal_acceleration = acceleration[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(solver::ExplicitSolver, model::SolidMechanics, integrator::CentralDifference)
    displacement = integrator.displacement
    velocity = integrator.velocity
    acceleration = solver.solution
    integrator.acceleration = acceleration
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = displacement[3*node-2:3*node]
        nodal_velocity = velocity[3*node-2:3*node]
        nodal_acceleration = acceleration[3*node-2:3*node]
        model.current[:, node] = model.reference[:, node] + nodal_displacement
        model.velocity[:, node] = nodal_velocity
        model.acceleration[:, node] = nodal_acceleration
    end
end

function copy_solution_source_targets(model::SolidMechanics, integrator::CentralDifference, solver::ExplicitSolver)
    _, num_nodes = size(model.reference)
    for node ∈ 1:num_nodes
        nodal_displacement = model.current[:, node] - model.reference[:, node]
        nodal_velocity = model.velocity[:, node]
        nodal_acceleration = model.acceleration[:, node]
        integrator.displacement[3*node-2:3*node] = nodal_displacement
        integrator.velocity[3*node-2:3*node] = nodal_velocity
        integrator.acceleration[3*node-2:3*node] = nodal_acceleration
    end
    solver.solution = integrator.acceleration
end

function evaluate(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    strain_energy, internal_force, body_force, stiffness_matrix = evaluate(integrator, model)
    solver.value = strain_energy
    external_force = body_force + model.boundary_tractions_force
    solver.gradient = internal_force - external_force
    solver.hessian = stiffness_matrix
end

function evaluate(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    strain_energy, internal_force, body_force, stiffness_matrix, mass_matrix = evaluate(integrator, model)
    β = integrator.β
    Δt = integrator.time_step
    inertial_force = mass_matrix * integrator.acceleration
    kinetic_energy = 0.5 * dot(integrator.velocity, mass_matrix, integrator.velocity)
    external_force = body_force + model.boundary_tractions_force
    solver.hessian = stiffness_matrix + mass_matrix / β / Δt / Δt
    solver.value = strain_energy - external_force ⋅ integrator.displacement + kinetic_energy
    solver.gradient = internal_force - external_force + inertial_force
end

function evaluate(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    strain_energy, internal_force, body_force, lumped_mass = evaluate(integrator, model)
    inertial_force = lumped_mass .* integrator.acceleration
    kinetic_energy = 0.5 * lumped_mass ⋅ (integrator.velocity .* integrator.velocity)
    external_force = body_force + model.boundary_tractions_force
    solver.value = strain_energy - external_force ⋅ integrator.displacement + kinetic_energy
    solver.gradient = internal_force - external_force + inertial_force
    solver.lumped_hessian = lumped_mass
end

function compute_step(solver::HessianMinimizer, _::NewtonStep, free::BitVector)
    step = zeros(length(solver.gradient))
    step[free] = -solver.hessian[free, free] \ solver.gradient[free]
    return step
end

function compute_step(solver::ExplicitSolver, _::ExplicitStep, free::BitVector)
    step = zeros(length(solver.gradient))
    step[free] = -solver.gradient[free] ./ solver.lumped_hessian[free]
    return step
end

function update_convergence_criterion(solver::HessianMinimizer, absolute_error::Float64)
    solver.absolute_error = absolute_error
    solver.relative_error = solver.initial_norm > 0.0 ? absolute_error / solver.initial_norm : 0.0
    converged_absolute = solver.absolute_error ≤ solver.absolute_tolerance
    converged_relative = solver.relative_error ≤ solver.relative_tolerance
    solver.converged = converged_absolute || converged_relative
end

function update_convergence_criterion(solver::ExplicitSolver, _::Float64)
    solver.converged = true
end

function continue_solve(solver::HessianMinimizer, iteration_number::Int64)
    if solver.failed == true
        return false
    end
    zero_residual = solver.absolute_error == 0.0
    if zero_residual == true
        return false
    end
    exceeds_minimum_iterations = iteration_number > solver.minimum_iterations
    if exceeds_minimum_iterations == false
        return true
    end
    exceeds_maximum_iterations = iteration_number > solver.maximum_iterations
    if exceeds_maximum_iterations == true
        return false
    end
    continue_solving = solver.converged == false
    return continue_solving
end

function continue_solve(_::ExplicitSolver, _::Int64)
    return false
end

function solve(integrator::Any, model::SolidMechanics, solver::Any)
    predict(integrator, solver, model)
    evaluate(integrator, solver, model)
    residual = solver.gradient
    norm_residual = norm(residual[model.free_dofs])
    solver.initial_norm = norm_residual
    iteration_number = 0
    solver.failed = solver.failed || model.failed
    step_type = solver.step
    while true
        step = compute_step(solver, step_type, model.free_dofs)
        solver.solution += step
        correct(integrator, solver, model)
        evaluate(integrator, solver, model)
        residual = solver.gradient
        norm_step = norm(step)
        norm_residual = norm(residual[model.free_dofs])
        if iteration_number == 0
            println("initial |R|=", norm_residual, ", |X|=", norm(solver.solution))
        else
            println("iter=", iteration_number, ", |R|=", norm_residual, ", |X|=", norm(solver.solution), ", |ΔX|=", norm_step)
        end
        update_convergence_criterion(solver, norm_residual)
        iteration_number += 1
        if continue_solve(solver, iteration_number) == false
            break
        end
    end
end