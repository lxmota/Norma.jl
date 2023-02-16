using DelimitedFiles
using Formatting

include("solver_def.jl")

function QuasiStatic(params::Dict{Any,Any})
    time_integrator_params = params["time integrator"]
    initial_time = time_integrator_params["initial time"]
    final_time = time_integrator_params["final time"]
    time_step = time_integrator_params["time step"]
    time = initial_time
    stop = 0
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    displacement = zeros(num_dof)
    QuasiStatic(initial_time, final_time, time_step, time, stop, displacement)
end

function Newmark(params::Dict{Any,Any})
    time_integrator_params = params["time integrator"]
    initial_time = time_integrator_params["initial time"]
    final_time = time_integrator_params["final time"]
    time_step = time_integrator_params["time step"]
    time = initial_time
    stop = 0
    β = time_integrator_params["β"]
    γ = time_integrator_params["γ"]
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    disp_pre = zeros(num_dof)
    velo_pre = zeros(num_dof)
    Newmark(initial_time, final_time, time_step, time, stop, β, γ, displacement, velocity, acceleration, disp_pre, velo_pre)
end

function CentralDifference(params::Dict{Any,Any})
    time_integrator_params = params["time integrator"]
    initial_time = time_integrator_params["initial time"]
    final_time = time_integrator_params["final time"]
    time_step = time_integrator_params["time step"]
    time = initial_time
    stop = 0
    γ = time_integrator_params["γ"]
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    CentralDifference(initial_time, final_time, time_step, time, stop, γ, displacement, velocity, acceleration)
end

function create_time_integrator(params::Dict{Any,Any})
    time_integrator_params = params["time integrator"]
    time_integrator_name = time_integrator_params["type"]
    if time_integrator_name == "quasi static"
        return QuasiStatic(params)
    elseif time_integrator_name == "Newmark"
        return Newmark(params)
    elseif time_integrator_name == "central difference"
        return CentralDifference(params)
    else
        error("Unknown type of time integrator : ", time_integrator_name)
    end
end

function predict(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    copy_solution_source_targets(model, integrator, solver)
end

function correct(integrator::QuasiStatic, solver::HessianMinimizer, model::SolidMechanics)
    copy_solution_source_targets(solver, model, integrator)
end

function predict(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    copy_solution_source_targets(model, integrator, solver)
    free = solver.free_dofs
    fixed = solver.free_dofs .== false
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    if integrator.stop == 0
        _, internal_force, external_force, _, mass_matrix = evaluate(integrator, model)
        inertial_force = external_force - internal_force
        integrator.acceleration[free] = mass_matrix[free, free] \ inertial_force[free]
    end
    u = integrator.displacement
    v = integrator.velocity
    a = integrator.acceleration
    uᵖʳᵉ = integrator.disp_pre
    vᵖʳᵉ = integrator.velo_pre
    uᵖʳᵉ[fixed] = u[fixed]
    vᵖʳᵉ[fixed] = v[fixed]
    uᵖʳᵉ[free] = u[free] + Δt * v[free] + (0.5 - β) * Δt * Δt * a[free]
    vᵖʳᵉ[free] = v[free] + (1.0 - γ) * Δt * a[free]
    if integrator.stop > 0
        u[free] = uᵖʳᵉ[free]
        v[free] = vᵖʳᵉ[free]
    end
    copy_solution_source_targets(integrator, solver, model)
end

function correct(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    u = integrator.displacement = solver.solution
    uᵖʳᵉ = integrator.disp_pre
    vᵖʳᵉ = integrator.velo_pre
    free = solver.free_dofs
    if integrator.stop > 0
        integrator.acceleration[free] = (u[free] - uᵖʳᵉ[free]) / β / Δt / Δt
        integrator.velocity[free] = vᵖʳᵉ[free] + γ * Δt * (u[free] - uᵖʳᵉ[free]) / β / Δt / Δt
    end
    copy_solution_source_targets(integrator, solver, model)
end

function predict(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    copy_solution_source_targets(model, integrator, solver)
    free = solver.free_dofs
    Δt = integrator.time_step
    γ = integrator.γ
    if integrator.stop == 0
        _, internal_force, external_force, lumped_mass = evaluate(integrator, model)
        inertial_force = external_force - internal_force
        integrator.acceleration[free] = inertial_force[free] ./ lumped_mass[free]
    else
        u = integrator.displacement
        v = integrator.velocity
        a = integrator.acceleration
        u[free] += Δt * v[free] + 0.5 * Δt * Δt * a[free]
        v[free] += (1.0 - γ) * Δt * a[free]
    end
    copy_solution_source_targets(integrator, solver, model)
end

function correct(integrator::CentralDifference, solver::ExplicitSolver, model::SolidMechanics)
    Δt = integrator.time_step
    γ = integrator.γ
    a = integrator.acceleration = solver.solution
    free = solver.free_dofs
    if integrator.stop > 0
        integrator.velocity[free] += γ * Δt * a[free]
    end
    copy_solution_source_targets(integrator, solver, model)
end

function initialize_writing(_::StaticTimeIntegrator, model::SolidMechanics)
    output_mesh = model.params["output_mesh"]
    num_node_vars = output_mesh.get_node_variable_number()
    disp_x_index = num_node_vars + 1
    disp_y_index = num_node_vars + 2
    disp_z_index = num_node_vars + 3
    num_node_vars += 3
    refe_x_index = num_node_vars + 1
    refe_y_index = num_node_vars + 2
    refe_z_index = num_node_vars + 3
    num_node_vars += 3
    output_mesh.set_node_variable_number(num_node_vars)
    output_mesh.put_node_variable_name("refe_x", refe_x_index)
    output_mesh.put_node_variable_name("refe_y", refe_y_index)
    output_mesh.put_node_variable_name("refe_z", refe_z_index)
    output_mesh.put_node_variable_name("disp_x", disp_x_index)
    output_mesh.put_node_variable_name("disp_y", disp_y_index)
    output_mesh.put_node_variable_name("disp_z", disp_z_index)
    num_element_vars = output_mesh.get_element_variable_number()
    elem_blk_ids = output_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    max_num_int_points = 0
    for blk_index ∈ 1:num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = output_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        max_num_int_points = max(max_num_int_points, num_points)
    end
    ip_var_index = num_element_vars
    num_element_vars += 6 * max_num_int_points
    output_mesh.set_element_variable_number(num_element_vars)
    for point ∈ 1:max_num_int_points
        stress_xx_index = ip_var_index + 1
        stress_yy_index = ip_var_index + 2
        stress_zz_index = ip_var_index + 3
        stress_yz_index = ip_var_index + 4
        stress_xz_index = ip_var_index + 5
        stress_xy_index = ip_var_index + 6
        ip_var_index += 6
        ip_str = sprintf1("_%d", point)
        output_mesh.put_element_variable_name("stress_xx" * ip_str, stress_xx_index)
        output_mesh.put_element_variable_name("stress_yy" * ip_str, stress_yy_index)
        output_mesh.put_element_variable_name("stress_zz" * ip_str, stress_zz_index)
        output_mesh.put_element_variable_name("stress_yz" * ip_str, stress_yz_index)
        output_mesh.put_element_variable_name("stress_xz" * ip_str, stress_xz_index)
        output_mesh.put_element_variable_name("stress_xy" * ip_str, stress_xy_index)
    end
end

function initialize_writing(_::DynamicTimeIntegrator, model::SolidMechanics)
    output_mesh = model.params["output_mesh"]
    num_node_vars = output_mesh.get_node_variable_number()
    disp_x_index = num_node_vars + 1
    disp_y_index = num_node_vars + 2
    disp_z_index = num_node_vars + 3
    num_node_vars += 3
    refe_x_index = num_node_vars + 1
    refe_y_index = num_node_vars + 2
    refe_z_index = num_node_vars + 3
    num_node_vars += 3
    velo_x_index = num_node_vars + 1
    velo_y_index = num_node_vars + 2
    velo_z_index = num_node_vars + 3
    num_node_vars += 3
    acce_x_index = num_node_vars + 1
    acce_y_index = num_node_vars + 2
    acce_z_index = num_node_vars + 3
    num_node_vars += 3
    output_mesh.set_node_variable_number(num_node_vars)
    output_mesh.put_node_variable_name("refe_x", refe_x_index)
    output_mesh.put_node_variable_name("refe_y", refe_y_index)
    output_mesh.put_node_variable_name("refe_z", refe_z_index)
    output_mesh.put_node_variable_name("disp_x", disp_x_index)
    output_mesh.put_node_variable_name("disp_y", disp_y_index)
    output_mesh.put_node_variable_name("disp_z", disp_z_index)
    output_mesh.put_node_variable_name("velo_x", velo_x_index)
    output_mesh.put_node_variable_name("velo_y", velo_y_index)
    output_mesh.put_node_variable_name("velo_z", velo_z_index)
    output_mesh.put_node_variable_name("acce_x", acce_x_index)
    output_mesh.put_node_variable_name("acce_y", acce_y_index)
    output_mesh.put_node_variable_name("acce_z", acce_z_index)
    num_element_vars = output_mesh.get_element_variable_number()
    elem_blk_ids = output_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    max_num_int_points = 0
    for blk_index ∈ 1:num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = output_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        max_num_int_points = max(max_num_int_points, num_points)
    end
    ip_var_index = num_element_vars
    num_element_vars += 6 * max_num_int_points
    output_mesh.set_element_variable_number(num_element_vars)
    for point ∈ 1:max_num_int_points
        stress_xx_index = ip_var_index + 1
        stress_yy_index = ip_var_index + 2
        stress_zz_index = ip_var_index + 3
        stress_yz_index = ip_var_index + 4
        stress_xz_index = ip_var_index + 5
        stress_xy_index = ip_var_index + 6
        ip_var_index += 6
        ip_str = sprintf1("_%d", point)
        output_mesh.put_element_variable_name("stress_xx" * ip_str, stress_xx_index)
        output_mesh.put_element_variable_name("stress_yy" * ip_str, stress_yy_index)
        output_mesh.put_element_variable_name("stress_zz" * ip_str, stress_zz_index)
        output_mesh.put_element_variable_name("stress_yz" * ip_str, stress_yz_index)
        output_mesh.put_element_variable_name("stress_xz" * ip_str, stress_xz_index)
        output_mesh.put_element_variable_name("stress_xy" * ip_str, stress_xy_index)
    end
end

function finalize_writing(model::Any)
    input_mesh = model.params["input_mesh"]
    input_mesh.close()
    output_mesh = model.params["output_mesh"]
    output_mesh.close()
end

function write_step(integrator::Any, model::Any)
    stop = integrator.stop
    exodus_interval = 1
    if haskey(model.params, "Exodus output interval") == true
        exodus_interval = model.params["Exodus output interval"]
    end
    if exodus_interval > 0 && stop % exodus_interval == 0
        write_step_exodus(integrator, model)
    end
    csv_interval = 0
    if haskey(model.params, "CSV output interval") == true
        csv_interval = model.params["CSV output interval"]
    end
    if csv_interval > 0 && stop % csv_interval == 0
        write_step_csv(integrator)
    end
end

function write_step_csv(integrator::StaticTimeIntegrator)
    stop = integrator.stop
    index_string = "-" * string(stop, pad=4)
    disp_filename = "disp" * index_string * ".csv"
    writedlm(disp_filename, integrator.displacement, '\n')
end

function write_step_csv(integrator::DynamicTimeIntegrator)
    stop = integrator.stop
    index_string = "-" * string(stop, pad=4)
    disp_filename = "disp" * index_string * ".csv"
    velo_filename = "velo" * index_string * ".csv"
    acce_filename = "acce" * index_string * ".csv"
    writedlm(disp_filename, integrator.displacement, '\n')
    writedlm(velo_filename, integrator.velocity, '\n')
    writedlm(acce_filename, integrator.acceleration, '\n')
end

function write_step_exodus(integrator::StaticTimeIntegrator, model::SolidMechanics)
    time = integrator.time
    stop = integrator.stop
    time_index = stop + 1
    output_mesh = model.params["output_mesh"]
    output_mesh.put_time(time_index, time)
    displacement = model.current - model.reference
    refe_x = model.current[1, :]
    refe_y = model.current[2, :]
    refe_z = model.current[3, :]
    output_mesh.put_variable_values("EX_NODAL", 1, "refe_x", time_index, refe_x)
    output_mesh.put_variable_values("EX_NODAL", 1, "refe_y", time_index, refe_y)
    output_mesh.put_variable_values("EX_NODAL", 1, "refe_z", time_index, refe_z)
    disp_x = displacement[1, :]
    disp_y = displacement[2, :]
    disp_z = displacement[3, :]
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_x", time_index, disp_x)
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_y", time_index, disp_y)
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_z", time_index, disp_z)
    stress = model.stress
    elem_blk_ids = output_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1:num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = output_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        blk_conn = output_mesh.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        block_stress = stress[blk_index]
        stress_xx = zeros(num_blk_elems, num_points)
        stress_yy = zeros(num_blk_elems, num_points)
        stress_zz = zeros(num_blk_elems, num_points)
        stress_yz = zeros(num_blk_elems, num_points)
        stress_xz = zeros(num_blk_elems, num_points)
        stress_xy = zeros(num_blk_elems, num_points)
        for blk_elem_index ∈ 1:num_blk_elems
            element_stress = block_stress[blk_elem_index]
            for point ∈ 1:num_points
                point_stress = element_stress[point]
                stress_xx[blk_elem_index, point] = point_stress[1]
                stress_yy[blk_elem_index, point] = point_stress[2]
                stress_zz[blk_elem_index, point] = point_stress[3]
                stress_yz[blk_elem_index, point] = point_stress[4]
                stress_xz[blk_elem_index, point] = point_stress[5]
                stress_xy[blk_elem_index, point] = point_stress[6]
            end
        end
        for point ∈ 1:num_points
            ip_str = sprintf1("_%d", point)
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xx" * ip_str, time_index, stress_xx[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_yy" * ip_str, time_index, stress_yy[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_zz" * ip_str, time_index, stress_zz[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_yz" * ip_str, time_index, stress_yz[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xz" * ip_str, time_index, stress_xz[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xy" * ip_str, time_index, stress_xy[:, point])
        end
    end
end

function write_step_exodus(integrator::DynamicTimeIntegrator, model::SolidMechanics)
    time = integrator.time
    stop = integrator.stop
    time_index = stop + 1
    output_mesh = model.params["output_mesh"]
    output_mesh.put_time(time_index, time)
    displacement = model.current - model.reference
    refe_x = model.current[1, :]
    refe_y = model.current[2, :]
    refe_z = model.current[3, :]
    output_mesh.put_variable_values("EX_NODAL", 1, "refe_x", time_index, refe_x)
    output_mesh.put_variable_values("EX_NODAL", 1, "refe_y", time_index, refe_y)
    output_mesh.put_variable_values("EX_NODAL", 1, "refe_z", time_index, refe_z)
    disp_x = displacement[1, :]
    disp_y = displacement[2, :]
    disp_z = displacement[3, :]
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_x", time_index, disp_x)
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_y", time_index, disp_y)
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_z", time_index, disp_z)
    velocity = model.velocity
    velo_x = velocity[1, :]
    velo_y = velocity[2, :]
    velo_z = velocity[3, :]
    output_mesh.put_variable_values("EX_NODAL", 1, "velo_x", time_index, velo_x)
    output_mesh.put_variable_values("EX_NODAL", 1, "velo_y", time_index, velo_y)
    output_mesh.put_variable_values("EX_NODAL", 1, "velo_z", time_index, velo_z)
    acceleration = model.acceleration
    acce_x = acceleration[1, :]
    acce_y = acceleration[2, :]
    acce_z = acceleration[3, :]
    output_mesh.put_variable_values("EX_NODAL", 1, "acce_x", time_index, acce_x)
    output_mesh.put_variable_values("EX_NODAL", 1, "acce_y", time_index, acce_y)
    output_mesh.put_variable_values("EX_NODAL", 1, "acce_z", time_index, acce_z)
    stress = model.stress
    elem_blk_ids = output_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1:num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = output_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        blk_conn = output_mesh.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        block_stress = stress[blk_index]
        stress_xx = zeros(num_blk_elems, num_points)
        stress_yy = zeros(num_blk_elems, num_points)
        stress_zz = zeros(num_blk_elems, num_points)
        stress_yz = zeros(num_blk_elems, num_points)
        stress_xz = zeros(num_blk_elems, num_points)
        stress_xy = zeros(num_blk_elems, num_points)
        for blk_elem_index ∈ 1:num_blk_elems
            element_stress = block_stress[blk_elem_index]
            for point ∈ 1:num_points
                point_stress = element_stress[point]
                stress_xx[blk_elem_index, point] = point_stress[1]
                stress_yy[blk_elem_index, point] = point_stress[2]
                stress_zz[blk_elem_index, point] = point_stress[3]
                stress_yz[blk_elem_index, point] = point_stress[4]
                stress_xz[blk_elem_index, point] = point_stress[5]
                stress_xy[blk_elem_index, point] = point_stress[6]
            end
        end
        for point ∈ 1:num_points
            ip_str = sprintf1("_%d", point)
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xx" * ip_str, time_index, stress_xx[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_yy" * ip_str, time_index, stress_yy[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_zz" * ip_str, time_index, stress_zz[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_yz" * ip_str, time_index, stress_yz[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xz" * ip_str, time_index, stress_xz[:, point])
            output_mesh.put_variable_values("EX_ELEM_BLOCK", blk_id, "stress_xy" * ip_str, time_index, stress_xy[:, point])
        end
    end
end