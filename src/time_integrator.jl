using DelimitedFiles
using Format

function adaptive_stepping_parameters(integrator_params::Dict{Any,Any})
    has_minimum = haskey(integrator_params, "minimum time step")
    has_decrease = haskey(integrator_params, "decrease factor")
    has_maximum = haskey(integrator_params, "maximum time step")
    has_increase = haskey(integrator_params, "increase factor")
    has_any = has_minimum || has_decrease || has_maximum || has_increase
    has_all = has_minimum && has_decrease && has_maximum && has_increase
    if (has_any == true && has_all == false)
        error("Adaptive time stepping requires 4 parameters: minimum and maximum time steps and decrease and increase factors")
    elseif (has_any == true && has_all == true)
        minimum_time_step = integrator_params["minimum time step"]
        decrease_factor = integrator_params["decrease factor"]
        maximum_time_step = integrator_params["maximum time step"]
        increase_factor = integrator_params["increase factor"]
    else
        minimum_time_step = maximum_time_step = integrator_params["time step"]
        decrease_factor = increase_factor = 1.0
    end
    return minimum_time_step, decrease_factor, maximum_time_step, increase_factor
end

function QuasiStatic(params::Dict{Any,Any})
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(integrator_params)
    time = prev_time = initial_time
    stop = 0
    input_mesh = params["input_mesh"]
    num_nodes = Exodus.num_nodes(input_mesh.init)
    num_dof = 3 * num_nodes
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    stored_energy = 0.0
    QuasiStatic(
        initial_time,
        final_time,
        time_step,
        minimum_time_step,
        decrease_factor,
        maximum_time_step,
        increase_factor,
        prev_time,
        time,
        stop,
        displacement,
        velocity,
        acceleration,
        stored_energy,
    )
end


## Have this for OpInf right now, idealy work in to current newmark
function NewmarkGeneral(params::Dict{Any,Any},model::LinearOpInfRom)
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(integrator_params)
    time = prev_time = initial_time
    stop = 0
    β = integrator_params["β"]
    γ = integrator_params["γ"]
    input_mesh = params["input_mesh"]
    input_mesh = params["input_mesh"]
    num_dof, = size(model.reduced_state)
    state = zeros(num_dof)
    state_dot = zeros(num_dof)
    state_ddot = zeros(num_dof)
    state_np1 = zeros(num_dof)
    state_np1_dot = zeros(num_dof)
    state_np1_ddot = zeros(num_dof)
    disp_pre = zeros(num_dof)
    velo_pre = zeros(num_dof)
    stored_energy = 0.0
    kinetic_energy = 0.0
    NewmarkGeneral(
        initial_time,
        final_time,
        time_step,
        minimum_time_step,
        decrease_factor,
        maximum_time_step,
        increase_factor,
        prev_time,
        time,
        stop,
        β,
        γ,
        state,
        state_dot,
        state_ddot,
        state_np1,
        state_np1_dot,
        state_np1_ddot,
        stored_energy,
        kinetic_energy,
    )
end


function Newmark(params::Dict{Any,Any},model::SolidMechanics)
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = integrator_params["time step"]
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(integrator_params)
    time = prev_time = initial_time
    stop = 0
    β = integrator_params["β"]
    γ = integrator_params["γ"]
    input_mesh = params["input_mesh"]
    input_mesh = params["input_mesh"]
    num_nodes = Exodus.num_nodes(input_mesh.init)
    num_dof = 3 * num_nodes
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    disp_pre = zeros(num_dof)
    velo_pre = zeros(num_dof)
    stored_energy = 0.0
    kinetic_energy = 0.0
    Newmark(
        initial_time,
        final_time,
        time_step,
        minimum_time_step,
        decrease_factor,
        maximum_time_step,
        increase_factor,
        prev_time,
        time,
        stop,
        β,
        γ,
        displacement,
        velocity,
        acceleration,
        disp_pre,
        velo_pre,
        stored_energy,
        kinetic_energy,
    )
end

function CentralDifference(params::Dict{Any,Any})
    integrator_params = params["time integrator"]
    initial_time = integrator_params["initial time"]
    final_time = integrator_params["final time"]
    time_step = 0.0
    minimum_time_step, decrease_factor, maximum_time_step, increase_factor = adaptive_stepping_parameters(integrator_params)
    stable_time_step = 0.0
    user_time_step = integrator_params["time step"]
    time = prev_time = initial_time
    stop = 0
    CFL = integrator_params["CFL"]
    γ = integrator_params["γ"]
    input_mesh = params["input_mesh"]
    input_mesh = params["input_mesh"]
    num_nodes = Exodus.num_nodes(input_mesh.init)
    num_dof = 3 * num_nodes
    displacement = zeros(num_dof)
    velocity = zeros(num_dof)
    acceleration = zeros(num_dof)
    stored_energy = 0.0
    kinetic_energy = 0.0
    CentralDifference(
        initial_time,
        final_time,
        time_step,
        minimum_time_step,
        decrease_factor,
        maximum_time_step,
        increase_factor,
        user_time_step,
        stable_time_step,
        prev_time,
        time,
        stop,
        CFL,
        γ,
        displacement,
        velocity,
        acceleration,
        stored_energy,
        kinetic_energy,
    )
end

function create_time_integrator(params::Dict{Any,Any},model::Any)
    integrator_params = params["time integrator"]
    integrator_name = integrator_params["type"]
    if integrator_name == "quasi static"
        return QuasiStatic(params)
    elseif integrator_name == "Newmark"
        return Newmark(params,model)
    elseif integrator_name == "NewmarkGeneral"
        return NewmarkGeneral(params,model)
    elseif integrator_name == "central difference"
        return CentralDifference(params)
    else
        error("Unknown type of time integrator : ", integrator_name)
    end
end

function is_static_or_dynamic(integrator_name::String)
    if integrator_name == "quasi static"
        return "static"
    elseif integrator_name == "Newmark"
        return "dynamic"
    elseif integrator_name == "central difference"
        return "dynamic"
    else
        error("Unknown type of time integrator : ", integrator_name)
    end
end

## Operator Inference 
function initialize(integrator::NewmarkGeneral, solver::Any, model::LinearOpInfRom)
    integrator.state_np1 = 1.0 .* model.reduced_state 
    integrator.state = 1.0 .* model.reduced_state 

end

function predict(integrator::NewmarkGeneral, solver::Any, model::LinearOpInfRom)
    # Move updated states
    integrator.state = 1.0.*integrator.state_np1
    integrator.state_dot = 1.0.*integrator.state_np1_dot
    integrator.state_ddot = 1.0.*integrator.state_np1_ddot
    solver.solution = 1.0 .* integrator.state
    #copy_solution_source_targets(model, integrator, solver)
end

function correct(integrator::NewmarkGeneral, solver::Any, model::LinearOpInfRom)
    model.reduced_state[:] = solver.solution[:]
end
###


function initialize(_::QuasiStatic, _::Any, _::SolidMechanics)
end

function predict(integrator::QuasiStatic, solver::Any, model::SolidMechanics)
    copy_solution_source_targets(model, integrator, solver)
end

function correct(integrator::QuasiStatic, solver::Any, model::SolidMechanics)
    copy_solution_source_targets(solver, model, integrator)
end

function initialize(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    println("Computing initial acceleration")
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    stored_energy, internal_force, external_force, _, mass_matrix =
        evaluate(integrator, model)
    inertial_force = external_force - internal_force
    kinetic_energy = 0.5 * dot(integrator.velocity, mass_matrix, integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    integrator.stored_energy = stored_energy
    integrator.acceleration[free] = mass_matrix[free, free] \ inertial_force[free]
    copy_solution_source_targets(integrator, solver, model)
end

function predict(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    fixed = .!free
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    u = integrator.displacement
    v = integrator.velocity
    a = integrator.acceleration
    uᵖʳᵉ = integrator.disp_pre
    vᵖʳᵉ = integrator.velo_pre
    uᵖʳᵉ[free] = u[free] += Δt * v[free] + (0.5 - β) * Δt * Δt * a[free]
    vᵖʳᵉ[free] = v[free] += (1.0 - γ) * Δt * a[free]
    uᵖʳᵉ[fixed] = u[fixed]
    vᵖʳᵉ[fixed] = v[fixed]
    copy_solution_source_targets(integrator, solver, model)
end

function correct(integrator::Newmark, solver::HessianMinimizer, model::SolidMechanics)
    free = model.free_dofs
    Δt = integrator.time_step
    β = integrator.β
    γ = integrator.γ
    u = integrator.displacement = solver.solution
    uᵖʳᵉ = integrator.disp_pre
    vᵖʳᵉ = integrator.velo_pre
    integrator.acceleration[free] = (u[free] - uᵖʳᵉ[free]) / β / Δt / Δt
    integrator.velocity[free] = vᵖʳᵉ[free] + γ * Δt * integrator.acceleration[free]
    copy_solution_source_targets(integrator, solver, model)
end

function initialize(
    integrator::CentralDifference,
    solver::ExplicitSolver,
    model::SolidMechanics,
)
    println("Computing initial acceleration")
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    set_time_step(integrator, model)
    stored_energy, internal_force, external_force, lumped_mass = evaluate(integrator, model)
    kinetic_energy = 0.5 * lumped_mass ⋅ (integrator.velocity .* integrator.velocity)
    integrator.kinetic_energy = kinetic_energy
    integrator.stored_energy = stored_energy
    inertial_force = external_force - internal_force
    integrator.acceleration[free] = inertial_force[free] ./ lumped_mass[free]
    copy_solution_source_targets(integrator, solver, model)
end

function predict(
    integrator::CentralDifference,
    solver::ExplicitSolver,
    model::SolidMechanics,
)
    copy_solution_source_targets(model, integrator, solver)
    free = model.free_dofs
    set_time_step(integrator, model)
    Δt = integrator.time_step
    γ = integrator.γ
    u = integrator.displacement
    v = integrator.velocity
    a = integrator.acceleration
    u[free] += Δt * v[free] + 0.5 * Δt * Δt * a[free]
    v[free] += (1.0 - γ) * Δt * a[free]
    copy_solution_source_targets(integrator, solver, model)
end

function correct(
    integrator::CentralDifference,
    solver::ExplicitSolver,
    model::SolidMechanics,
)
    Δt = integrator.time_step
    γ = integrator.γ
    a = integrator.acceleration = solver.solution
    free = model.free_dofs
    integrator.velocity[free] += γ * Δt * a[free]
    copy_solution_source_targets(integrator, solver, model)
end



function initialize_writing(
    params::Dict{Any,Any},
    _::StaticTimeIntegrator,
    model::HeatConduction,
)
    output_mesh = params["output_mesh"]
    num_node_vars = Exodus.read_number_of_variables(output_mesh, NodalVariable)
    temp_index = num_node_vars + 1
    num_node_vars += 1
    Exodus.write_number_of_variables(output_mesh, NodalVariable, num_node_vars)
    Exodus.write_name(output_mesh, NodalVariable, Int32(temp_index), "temp")
    num_element_vars = Exodus.read_number_of_variables(output_mesh, ElementVariable)
    blocks = Exodus.read_sets(output_mesh, Block)
    max_num_int_points = 0
    for block ∈ blocks
        blk_id = block.id
        element_type = Exodus.read_block_parameters(output_mesh, blk_id)[1]
        num_points = default_num_int_pts(element_type)
        max_num_int_points = max(max_num_int_points, num_points)
    end
    ip_var_index = num_element_vars
end

function initialize_writing(
    params::Dict{Any,Any},
    integrator::DynamicTimeIntegrator,
    model::LinearOpInfRom,
)
    initialize_writing(params,integrator,model.fom_model)
end


function initialize_writing(
    params::Dict{Any,Any},
    _::StaticTimeIntegrator,
    model::SolidMechanics,
)
    output_mesh = params["output_mesh"]
    num_node_vars = Exodus.read_number_of_variables(output_mesh, NodalVariable)
    disp_x_index = num_node_vars + 1
    disp_y_index = num_node_vars + 2
    disp_z_index = num_node_vars + 3
    num_node_vars += 3
    refe_x_index = num_node_vars + 1
    refe_y_index = num_node_vars + 2
    refe_z_index = num_node_vars + 3
    num_node_vars += 3
    Exodus.write_number_of_variables(output_mesh, NodalVariable, num_node_vars)
    Exodus.write_name(output_mesh, NodalVariable, Int32(refe_x_index), "refe_x")
    Exodus.write_name(output_mesh, NodalVariable, Int32(refe_y_index), "refe_y")
    Exodus.write_name(output_mesh, NodalVariable, Int32(refe_z_index), "refe_z")
    Exodus.write_name(output_mesh, NodalVariable, Int32(disp_x_index), "disp_x")
    Exodus.write_name(output_mesh, NodalVariable, Int32(disp_y_index), "disp_y")
    Exodus.write_name(output_mesh, NodalVariable, Int32(disp_z_index), "disp_z")
    num_element_vars = Exodus.read_number_of_variables(output_mesh, ElementVariable)
    blocks = Exodus.read_sets(output_mesh, Block)
    max_num_int_points = 0
    for block ∈ blocks
        blk_id = block.id
        element_type = Exodus.read_block_parameters(output_mesh, blk_id)[1]
        num_points = default_num_int_pts(element_type)
        max_num_int_points = max(max_num_int_points, num_points)
    end
    ip_var_index = num_element_vars
    num_element_vars += 6 * max_num_int_points
    num_element_vars += 1
    Exodus.write_number_of_variables(output_mesh, ElementVariable, num_element_vars)
    for point ∈ 1:max_num_int_points
        stress_xx_index = ip_var_index + 1
        stress_yy_index = ip_var_index + 2
        stress_zz_index = ip_var_index + 3
        stress_yz_index = ip_var_index + 4
        stress_xz_index = ip_var_index + 5
        stress_xy_index = ip_var_index + 6
        ip_var_index += 6
        ip_str = "_" * cfmt("%d", point)
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_xx_index),
            "stress_xx" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_yy_index),
            "stress_yy" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_zz_index),
            "stress_zz" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_yz_index),
            "stress_yz" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_xz_index),
            "stress_xz" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_xy_index),
            "stress_xy" * ip_str,
        )
    end
    element_stored_energy_index = ip_var_index
    Exodus.write_name(
        output_mesh,
        ElementVariable,
        Int32(element_stored_energy_index),
        "stored_energy",
    )
end

function initialize_writing(
    params::Dict{Any,Any},
    _::DynamicTimeIntegrator,
    model::SolidMechanics,
)
    output_mesh = params["output_mesh"]
    num_node_vars = Exodus.read_number_of_variables(output_mesh, NodalVariable)
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
    Exodus.write_number_of_variables(output_mesh, NodalVariable, num_node_vars)
    Exodus.write_name(output_mesh, NodalVariable, Int32(refe_x_index), "refe_x")
    Exodus.write_name(output_mesh, NodalVariable, Int32(refe_y_index), "refe_y")
    Exodus.write_name(output_mesh, NodalVariable, Int32(refe_z_index), "refe_z")
    Exodus.write_name(output_mesh, NodalVariable, Int32(disp_x_index), "disp_x")
    Exodus.write_name(output_mesh, NodalVariable, Int32(disp_y_index), "disp_y")
    Exodus.write_name(output_mesh, NodalVariable, Int32(disp_z_index), "disp_z")
    Exodus.write_name(output_mesh, NodalVariable, Int32(velo_x_index), "velo_x")
    Exodus.write_name(output_mesh, NodalVariable, Int32(velo_y_index), "velo_y")
    Exodus.write_name(output_mesh, NodalVariable, Int32(velo_z_index), "velo_z")
    Exodus.write_name(output_mesh, NodalVariable, Int32(acce_x_index), "acce_x")
    Exodus.write_name(output_mesh, NodalVariable, Int32(acce_y_index), "acce_y")
    Exodus.write_name(output_mesh, NodalVariable, Int32(acce_z_index), "acce_z")
    num_element_vars = Exodus.read_number_of_variables(output_mesh, ElementVariable)
    blocks = Exodus.read_sets(output_mesh, Block)
    max_num_int_points = 0
    for block ∈ blocks
        blk_id = block.id
        element_type = Exodus.read_block_parameters(output_mesh, blk_id)[1]
        num_points = default_num_int_pts(element_type)
        max_num_int_points = max(max_num_int_points, num_points)
    end
    ip_var_index = num_element_vars
    num_element_vars += 6 * max_num_int_points
    num_element_vars += 1
    Exodus.write_number_of_variables(output_mesh, ElementVariable, num_element_vars)
    for point ∈ 1:max_num_int_points
        stress_xx_index = ip_var_index + 1
        stress_yy_index = ip_var_index + 2
        stress_zz_index = ip_var_index + 3
        stress_yz_index = ip_var_index + 4
        stress_xz_index = ip_var_index + 5
        stress_xy_index = ip_var_index + 6
        ip_var_index += 6
        ip_str = "_" * cfmt("%d", point)
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_xx_index),
            "stress_xx" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_yy_index),
            "stress_yy" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_zz_index),
            "stress_zz" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_yz_index),
            "stress_yz" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_xz_index),
            "stress_xz" * ip_str,
        )
        Exodus.write_name(
            output_mesh,
            ElementVariable,
            Int32(stress_xy_index),
            "stress_xy" * ip_str,
        )
    end
    element_stored_energy_index = ip_var_index
    Exodus.write_name(
        output_mesh,
        ElementVariable,
        Int32(element_stored_energy_index),
        "stored_energy",
    )
end

function finalize_writing(params::Dict{Any,Any})
    input_mesh = params["input_mesh"]
    Exodus.close(input_mesh)
    output_mesh = params["output_mesh"]
    Exodus.close(output_mesh)
end


function writedlm_nodal_array(filename::String, nodal_array::Matrix{Float64})
    open(filename, "w") do io
        for col ∈ 1:size(nodal_array, 2)
            # Write each column as a comma-separated line
            println(io, join(nodal_array[:, col], ","))
        end
    end
end

function write_step(params::Dict{Any,Any}, integrator::Any, model::Any)
    stop = integrator.stop
    exodus_interval = get(params, "Exodus output interval", 1)
    if exodus_interval > 0 && stop % exodus_interval == 0
        write_step_exodus(params, integrator, model)
    end
    csv_interval = get(params, "CSV output interval", 0)
    if csv_interval > 0 && stop % csv_interval == 0
        sim_id = 1
        if haskey(params, "global_simulation") == true
            sim_id = params["global_simulation"].subsim_name_index_map[params["name"]]
        end
        write_step_csv(integrator, model, sim_id)
        if haskey(params, "CSV write sidesets") == true
          write_sideset_step_csv(params,integrator,model,sim_id)
        end
    end
end

function write_step_csv(integrator::StaticTimeIntegrator, model::SolidMechanics, sim_id::Integer)
    stop = integrator.stop
    index_string = "-" * string(stop, pad = 4)
    sim_id_string = string(sim_id, pad = 2) * "-"
    curr_filename = sim_id_string * "curr" * index_string * ".csv"
    disp_filename = sim_id_string * "disp" * index_string * ".csv"
    potential_filename = sim_id_string * "potential" * index_string * ".csv"
    time_filename = sim_id_string * "time" * index_string * ".csv"
    writedlm_nodal_array(curr_filename, model.current)
    writedlm_nodal_array(disp_filename, model.current - model.reference)
    writedlm(potential_filename, integrator.stored_energy, '\n')
    writedlm(time_filename, integrator.time, '\n')
    if stop == 0
        refe_filename = sim_id_string * "refe" * ".csv"
        writedlm_nodal_array(refe_filename, model.reference)
    end
end


function write_sideset_step_csv(params::Dict{Any,Any},integrator::DynamicTimeIntegrator, model::SolidMechanics, sim_id::Integer)
  stop = integrator.stop
  index_string = "-" * string(stop, pad = 4)
  sim_id_string = string(sim_id, pad = 2) * "-"
  input_mesh = params["input_mesh"]
  bc_params = params["boundary conditions"]

  for bc ∈ model.boundary_conditions
      node_set_name = bc.node_set_name 
      curr_filename = sim_id_string * node_set_name * "curr" * index_string * ".csv"
      disp_filename = sim_id_string * node_set_name * "disp" * index_string * ".csv"
      velo_filename = sim_id_string * node_set_name * "velo" * index_string * ".csv"
      acce_filename = sim_id_string * node_set_name * "acce" * index_string * ".csv"

      writedlm(curr_filename, model.current[bc.offset,bc.node_set_node_indices])
      writedlm(velo_filename, model.velocity[bc.offset,bc.node_set_node_indices])
      writedlm(acce_filename, model.acceleration[bc.offset,bc.node_set_node_indices])
      writedlm(disp_filename, model.current[bc.offset,bc.node_set_node_indices] - model.reference[bc.offset,bc.node_set_node_indices])

  end



  #=
  for (bc_type, bc_type_params) ∈ bc_params
    for bc_setting_params ∈ bc_type_params
      node_set_name = bc_setting_params["node set"]
      node_set_id = node_set_id_from_name(node_set_name, input_mesh)
      side_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)
      curr_filename = sim_id_string * node_set_name * "curr" * index_string * ".csv"
      disp_filename = sim_id_string * node_set_name * "disp" * index_string * ".csv"
      velo_filename = sim_id_string * node_set_name * "velo" * index_string * ".csv"
      acce_filename = sim_id_string * node_set_name * "acce" * index_string * ".csv"
      writedlm(curr_filename, model.current[bc.offset,bc.node_set_node_indices])
      writedlm(velo_filename, model.velocity[bc.offset,side_set_node_indices])
      writedlm(acce_filename, model.acceleration[bc.offset,side_set_node_indices])
      writedlm(disp_filename, model.current[bc.offset,side_set_node_indices] - model.reference[bc.offset,side_set_node_indices])
    end
  end
  =#
end


function write_step_csv(integrator::DynamicTimeIntegrator, model::SolidMechanics, sim_id::Integer)
    stop = integrator.stop
    index_string = "-" * string(stop, pad = 4)
    sim_id_string = string(sim_id, pad = 2) * "-"
    free_dofs_filename = sim_id_string * "free_dofs" * index_string * ".csv"
    curr_filename = sim_id_string * "curr" * index_string * ".csv"
    disp_filename = sim_id_string * "disp" * index_string * ".csv"
    velo_filename = sim_id_string * "velo" * index_string * ".csv"
    acce_filename = sim_id_string * "acce" * index_string * ".csv"
    time_filename = sim_id_string * "time" * index_string * ".csv"
    potential_filename = sim_id_string * "potential" * index_string * ".csv"
    kinetic_filename = sim_id_string * "kinetic" * index_string * ".csv"
    writedlm(free_dofs_filename, model.free_dofs)
    writedlm_nodal_array(curr_filename, model.current)
    writedlm_nodal_array(velo_filename, model.velocity)
    writedlm_nodal_array(acce_filename, model.acceleration)
    writedlm_nodal_array(disp_filename, model.current - model.reference)
    writedlm(potential_filename, integrator.stored_energy, '\n')
    writedlm(kinetic_filename, integrator.kinetic_energy, '\n')
    writedlm(time_filename, integrator.time, '\n')
    if stop == 0 
        refe_filename = sim_id_string * "refe" * ".csv"
        writedlm_nodal_array(refe_filename, model.reference)
    end 
end



function write_step_csv(integrator::NewmarkGeneral, model::OpInfModel, sim_id::Integer)
    stop = integrator.stop
    index_string = "-" * string(stop, pad = 4)
    sim_id_string = string(sim_id, pad = 2) * "-"
    reduced_states_filename = sim_id_string * "reduced_states" * index_string * ".csv"
    time_filename = sim_id_string * "time" * index_string * ".csv"
    writedlm(reduced_states_filename, model.reduced_state)
    write_step_csv(integrator,model.fom_model,sim_id)
end



function write_step_exodus(
    params::Dict{Any,Any},
    integrator::StaticTimeIntegrator,
    model::HeatConduction,
)
    time = integrator.time
    stop = integrator.stop
    time_index = stop + 1
    output_mesh = params["output_mesh"]
    Exodus.write_time(output_mesh, time_index, time)
    temp = model.temperature[:]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "temp", temp)
end


function write_step_exodus(
    params::Dict{Any,Any},
    integrator::DynamicTimeIntegrator,
    model::LinearOpInfRom,
)
    #Re-construct full state
    reduced_state = model.reduced_state[:]
    for i = 1 : size(model.fom_model.current)[2]
      x_dof_index = 3 * (i - 1) + 1 
      y_dof_index = 3 * (i - 1) + 2 
      z_dof_index = 3 * (i - 1) + 3 
      if model.fom_model.free_dofs[x_dof_index]
        model.fom_model.current[1,i] = model.basis[1,i,:]'reduced_state + model.fom_model.reference[1,i]
      end  

      if model.fom_model.free_dofs[y_dof_index]
        model.fom_model.current[2,i] = model.basis[2,i,:]'reduced_state + model.fom_model.reference[2,i]
      end
 
      if model.fom_model.free_dofs[z_dof_index]
        model.fom_model.current[3,i] = model.basis[3,i,:]'reduced_state+ model.fom_model.reference[3,i]
      end
    end
    write_step_exodus(params,integrator,model.fom_model)
end


function write_step_exodus(
    params::Dict{Any,Any},
    integrator::StaticTimeIntegrator,
    model::SolidMechanics,
)
    time = integrator.time
    stop = integrator.stop
    time_index = stop + 1
    output_mesh = params["output_mesh"]
    Exodus.write_time(output_mesh, time_index, time)
    displacement = model.current - model.reference
    refe_x = model.current[1, :]
    refe_y = model.current[2, :]
    refe_z = model.current[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_x", refe_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_y", refe_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_z", refe_z)
    disp_x = displacement[1, :]
    disp_y = displacement[2, :]
    disp_z = displacement[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_x", disp_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_y", disp_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_z", disp_z)
    stress = model.stress
    stored_energy = model.stored_energy
    blocks = Exodus.read_sets(output_mesh, Block)
    for (block, block_stress, block_stored_energy) ∈ zip(blocks, stress, stored_energy)
        blk_id = block.id
        element_type, num_blk_elems, _, _, _, _ =
            Exodus.read_block_parameters(output_mesh, blk_id)
        num_points = default_num_int_pts(element_type)
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
            ip_str = "_" * cfmt("%d", point)
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_xx" * ip_str,
                stress_xx[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_yy" * ip_str,
                stress_yy[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_zz" * ip_str,
                stress_zz[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_yz" * ip_str,
                stress_yz[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_xz" * ip_str,
                stress_xz[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_xy" * ip_str,
                stress_xy[:, point],
            )
        end
        Exodus.write_values(
            output_mesh,
            ElementVariable,
            time_index,
            Int64(blk_id),
            "stored_energy",
            block_stored_energy,
        )
    end
end

function write_step_exodus(
    params::Dict{Any,Any},
    integrator::DynamicTimeIntegrator,
    model::SolidMechanics,
)
    time = integrator.time
    stop = integrator.stop
    time_index = stop + 1
    output_mesh = params["output_mesh"]
    Exodus.write_time(output_mesh, time_index, time)
    displacement = model.current - model.reference
    refe_x = model.current[1, :]
    refe_y = model.current[2, :]
    refe_z = model.current[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_x", refe_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_y", refe_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "refe_z", refe_z)
    disp_x = displacement[1, :]
    disp_y = displacement[2, :]
    disp_z = displacement[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_x", disp_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_y", disp_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "disp_z", disp_z)
    velocity = model.velocity
    velo_x = velocity[1, :]
    velo_y = velocity[2, :]
    velo_z = velocity[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_x", velo_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_y", velo_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "velo_z", velo_z)
    acceleration = model.acceleration
    acce_x = acceleration[1, :]
    acce_y = acceleration[2, :]
    acce_z = acceleration[3, :]
    Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_x", acce_x)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_y", acce_y)
    Exodus.write_values(output_mesh, NodalVariable, time_index, "acce_z", acce_z)
    stress = model.stress
    blocks = Exodus.read_sets(output_mesh, Block)
    for (block, block_stress) ∈ zip(blocks, stress)
        blk_id = block.id
        element_type, num_blk_elems, _, _, _, _ =
            Exodus.read_block_parameters(output_mesh, blk_id)
        num_points = default_num_int_pts(element_type)
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
            ip_str = "_" * cfmt("%d", point)
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_xx" * ip_str,
                stress_xx[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_yy" * ip_str,
                stress_yy[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_zz" * ip_str,
                stress_zz[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_yz" * ip_str,
                stress_yz[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_xz" * ip_str,
                stress_xz[:, point],
            )
            Exodus.write_values(
                output_mesh,
                ElementVariable,
                time_index,
                Int64(blk_id),
                "stress_xy" * ip_str,
                stress_xy[:, point],
            )
        end
    end
end
