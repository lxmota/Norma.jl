using Formatting
using SparseArrays

include("constitutive.jl")

abstract type Model end

@enum DOF free Dirichlet Schwarz

mutable struct SolidMechanics <: Model
    params::Dict{Any, Any}
    materials::Vector{Solid}
    reference::Matrix{Float64}
    current::Matrix{Float64}
    velocity::Matrix{Float64}
    acceleration::Matrix{Float64}
    stress::Vector{Vector{Vector{Vector{Float64}}}}
    nodal_dofs::Vector{DOF}
    time::Float64
    failed::Bool
end

function SolidMechanics(params::Dict{Any, Any})
    input_mesh = params["input_mesh"]
    x, y, z = input_mesh.get_coords()
    num_nodes = length(x)
    reference = Matrix{Float64}(undef, 3, num_nodes)
    current = Matrix{Float64}(undef, 3, num_nodes)
    velocity = Matrix{Float64}(undef, 3, num_nodes)
    acceleration = Matrix{Float64}(undef, 3, num_nodes)
    for node ∈ 1 : num_nodes
        reference[:, node] = [x[node], y[node], z[node]]
        current[:, node] = [x[node], y[node], z[node]]
        velocity[:, node] = [0.0, 0.0, 0.0]
        acceleration[:, node] = [0.0, 0.0, 0.0]
    end
    materials_file = params["material"]
    material_params = YAML.load_file(materials_file)
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    if (num_blks_params ≠ num_blks)
        error("number of blocks in mesh ", params["mesh"], " (", num_blks,
        ") must be equal to number of blocks in materials file ", params["material"],
        " (", num_blks_params, ")")
    end
    elem_blk_names = input_mesh.get_elem_blk_names()
    materials = Vector{Solid}(undef, 0)
    for elem_blk_name ∈ elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        push!(materials, material_model)
    end
    time = 0.0
    failed = false
    nodal_dofs = [free::DOF for _ ∈ 1 : 3 * num_nodes]
    stress = Vector{Vector{Vector{Vector{Float64}}}}(undef, num_blks)
    for blk_index ∈ 1 : num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        blk_conn = input_mesh.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        block_stress = Vector{Vector{Vector{Float64}}}(undef, num_blk_elems)
        for blk_elem_index ∈ 1 : num_blk_elems
            element_stress = Vector{Vector{Float64}}(undef, num_points)
            for point ∈ 1 : num_points
                element_stress[point] = zeros(6)
            end
            block_stress[blk_elem_index] = element_stress
        end
        stress[blk_index] = block_stress
    end
    SolidMechanics(params, materials, reference, current, velocity, acceleration, stress, nodal_dofs, time, failed)
end

mutable struct HeatConduction <: Model
    params::Dict{Any, Any}
    materials::Vector{Vector}
    reference::Matrix{Float64}
    temperature::Vector{Float64}
    rate::Vector{Float64}
    flux::Vector{Vector{Vector{Vector{Float64}}}}
    nodal_dofs::Vector{DOF}
    time::Float64
    failed::Bool
end

function HeatConduction(params::Dict{Any, Any})
    input_mesh = params["input_mesh"]
    x, y, z = input_mesh.get_coords()
    num_nodes = length(x)
    reference = Matrix{Float64}(undef, 3, num_nodes)
    temperature = Vector{Float64}(undef, num_nodes)
    rate = Vector{Float64}(undef, num_nodes)
    for node ∈ 1 : num_nodes
        reference[:, node] = [x[node], y[node], z[node]]
        temperature[node] = 0.0
        rate[node] = 0.0
    end
    materials_file = params["material"]
    material_params = YAML.load_file(materials_file)
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    if (num_blks_params ≠ num_blks)
        error("number of blocks in mesh ", params["mesh"], " (", num_blks,
        ") must be equal to number of blocks in materials file ", params["material"],
        " (", num_blks_params, ")")
    end
    elem_blk_names = input_mesh.get_elem_blk_names()
    materials = Vector{Thermal}(undef, 0)
    for elem_blk_name ∈ elem_blk_names
        material_name = material_blocks[elem_blk_name]
        material_props = material_params[material_name]
        material_model = create_material(material_props)
        push!(materials, material_model)
    end
    time = 0.0
    failed = false
    nodal_dofs = [free::DOF for _ ∈ 1 : num_nodes]
    flux = Vector{Vector{Vector{Vector{Float64}}}}(undef, num_blks)
    for blk_index ∈ 1 : num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        blk_conn = input_mesh.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        block_flux = Vector{Vector{Vector{Float64}}}(undef, num_blk_elems)
        for blk_elem_index ∈ 1 : num_blk_elems
            element_flux = Vector{Vector{Float64}}(undef, num_points)
            for point ∈ 1 : num_points
                element_flux[point] = zeros(3)
            end
            block_flux[blk_elem_index] = element_stress
        end
        flux[blk_index] = block_flux
    end
    HeatConduction(params, materials, reference, temperature, rate, flux, nodal_dofs, time, failed)
end

function create_model(params::Dict{Any, Any})
    model_name = params["model"]
    if model_name == "solid mechanics"
        return SolidMechanics(params)
    elseif model_name == "heat conduction"
        return HeatConduction(params)
    else
        error("Unknown type of model : ", model_name)
    end
end

function voigt_cauchy_from_stress(material::Any, P::MTTensor, F::MTTensor, J::Float64)
    σ = F * P' ./ J
    return [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]
end

function voigt_cauchy_from_stress(material::Linear_Elastic, σ::MTTensor, F::MTTensor, J::Float64)
    return [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]
end

function energy_force_stiffness(model::SolidMechanics)
    params = model.params
    materials = model.materials
    input_mesh = params["input_mesh"]
    solver_struct = params["solver_struct"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    total_energy = solver_struct.value = 0.0
    total_internal_force = solver_struct.gradient
    total_stiffness = solver_struct.hessian
    fill!(total_internal_force, 0.0)
    fill!(total_stiffness, 0.0)
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1 : num_blks
        material = materials[blk_index]
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        _, dNdξ, elem_weights = isoparametric(elem_type, num_points)
        elem_blk_conn, num_blk_elems, num_elem_nodes = input_mesh.get_elem_connectivity(blk_id)
        elem_dofs = zeros(Int64, 3 * num_elem_nodes)
        for blk_elem_index ∈ 1 : num_blk_elems
            conn_indices = (blk_elem_index - 1) * num_elem_nodes + 1 : blk_elem_index * num_elem_nodes
            node_indices = elem_blk_conn[conn_indices]
            elem_ref_pos = model.reference[:, node_indices]
            elem_cur_pos = model.current[:, node_indices]
            element_energy = 0.0
            element_internal_force = zeros(3 * num_elem_nodes)
            element_stiffness = zeros(3 * num_elem_nodes, 3 * num_elem_nodes)
            elem_dofs[1 : 3 : 3 * num_elem_nodes - 2] = 3 .* node_indices .- 2
            elem_dofs[2 : 3 : 3 * num_elem_nodes - 1] = 3 .* node_indices .- 1
            elem_dofs[3 : 3 : 3 * num_elem_nodes] = 3 .* node_indices
            for point ∈ 1 : num_points
                dXdξ = dNdξ[:, :, point] * elem_ref_pos'
                dxdξ = dNdξ[:, :, point] * elem_cur_pos'
                dxdX = dXdξ \ dxdξ
                dNdX = dXdξ \ dNdξ[:, :, point]
                B = gradient_operator(dNdX)
                j = det(dXdξ)
                J = det(dxdX)
                if J ≤ 0.0
                    model.failed = true
                    return 0.0, zeros(num_dof), spzeros(num_dof, num_dof)
                end
                F = MTTensor(dxdX)
                W, P, A = constitutive(material, F)
                stress = P[1:9]
                moduli = second_from_fourth(A)
                w = elem_weights[point]
                element_energy += W * j * w
                element_internal_force += B' * stress * j * w
                element_stiffness += B' * moduli * B * j * w
                voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
                model.stress[blk_index][blk_elem_index][point] = voigt_cauchy
            end
            total_energy += element_energy
            total_internal_force[elem_dofs] += element_internal_force
            total_stiffness[elem_dofs, elem_dofs] += element_stiffness
        end
    end
    return total_energy, total_internal_force, total_stiffness
end

function node_set_id_from_name(node_set_name::String, mesh::PyObject)
    node_set_names = mesh.get_node_set_names()
    num_names = length(node_set_names)
    node_set_index = 0
    for index ∈ 1 : num_names
        if (node_set_name == node_set_names[index])
            node_set_index = index
            break
        end
    end
    if (node_set_index == 0)
        error("node set ", node_set_name, " cannot be found in mesh")
    end
    node_set_ids = mesh.get_node_set_ids()
    node_set_id = node_set_ids[node_set_index]
    return node_set_id
end

function component_offset_from_string(name::String)
    offset = 0
    if name == "x"
        offset = 1
    elseif name == "y"
        offset = 2
    elseif name == "z"
        offset = 3
    else
        error("invalid component name ", name)
    end
    return offset
end

function apply_bcs(model::SolidMechanics)
    params = model.params
    reference = model.reference
    current = model.current
    input_mesh = params["input_mesh"]
    global t = model.time
    xc, yc, zc = input_mesh.get_coords()
    num_nodes = length(xc)
    nodal_dofs = [free::DOF for _ ∈ 1 : 3 * num_nodes]
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) ∈ bc_params
        for bc ∈ bc_type_params
            if bc_type == "Dirichlet"
                node_set_name = bc["node set"]
                function_str = bc["function"]
                component = bc["component"]
                offset = component_offset_from_string(component)
                node_set_id = node_set_id_from_name(node_set_name, input_mesh)
                node_set_node_indices = input_mesh.get_node_set_nodes(node_set_id)
                for node_index ∈ node_set_node_indices
                    global x = xc[node_index]
                    global y = yc[node_index]
                    global z = zc[node_index]
                    # function_str is an arbitrary function of t, x, y, z in the input file
                    bc_expr = Meta.parse(function_str)
                    bc_val = eval(bc_expr)
                    dof_index = 3 * (node_index - 1) + offset
                    current[dof_index] = reference[dof_index] + bc_val
                    nodal_dofs[dof_index] = Dirichlet::DOF
                end
            elseif bc_type == "Schwarz"
            elseif bc_type == "Neumann"
            end
        end
    end
    model.nodal_dofs = nodal_dofs
end

function apply_bcs(model::HeatConduction)
    params = model.params
    temperature = model.temperature
    input_mesh = params["input_mesh"]
    global t = model.time
    xc, yc, zc = input_mesh.get_coords()
    num_nodes = length(xc)
    nodal_dofs = [free::DOF for _ ∈ 1 : num_nodes]
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) ∈ bc_params
        for bc ∈ bc_type_params
            if bc_type == "Dirichlet"
                node_set_name = bc["node set"]
                function_str = bc["function"]
                node_set_id = node_set_id_from_name(node_set_name, input_mesh)
                node_set_node_indices = input_mesh.get_node_set_nodes(node_set_id)
                for node_index ∈ node_set_node_indices
                    global x = xc[node_index]
                    global y = yc[node_index]
                    global z = zc[node_index]
                    # function_str is an arbitrary function of t, x, y, z in the input file
                    bc_expr = Meta.parse(function_str)
                    bc_val = eval(bc_expr)
                    temperature[node_index] = bc_val
                    nodal_dofs[node_index] = Dirichlet::DOF
                end
            elseif bc_type == "Schwarz"
            elseif bc_type == "Neumann"
            end
        end
    end
    model.nodal_dofs = nodal_dofs
end

function initialize_writing(model::SolidMechanics)
    output_mesh = model.params["output_mesh"]
    num_node_vars = output_mesh.get_node_variable_number()
    disp_x_index = num_node_vars + 1
    disp_y_index = num_node_vars + 2
    disp_z_index = num_node_vars + 3
    num_node_vars += 3
    output_mesh.set_node_variable_number(num_node_vars)
    output_mesh.put_node_variable_name("disp_x", disp_x_index)
    output_mesh.put_node_variable_name("disp_y", disp_y_index)
    output_mesh.put_node_variable_name("disp_z", disp_z_index)
    num_element_vars = output_mesh.get_element_variable_number()
    elem_blk_ids = output_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    max_num_int_points = 0
    for blk_index ∈ 1 : num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = output_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        max_num_int_points = max(max_num_int_points, num_points)
    end
    ip_var_index = num_element_vars
    num_element_vars += 6 * max_num_int_points
    output_mesh.set_element_variable_number(num_element_vars)
    for point ∈ 1 : max_num_int_points
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
    output_mesh = model.params["output_mesh"]
    output_mesh.close()
end

function write_step(model::SolidMechanics, time_index::Int64, time::Float64)
    output_mesh = model.params["output_mesh"]
    output_mesh.put_time(time_index, time)
    displacement = model.current - model.reference
    disp_x = displacement[1, :]
    disp_y = displacement[2, :]
    disp_z = displacement[3, :]
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_x", time_index, disp_x)
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_y", time_index, disp_y)
    output_mesh.put_variable_values("EX_NODAL", 1, "disp_z", time_index, disp_z)
    stress = model.stress
    elem_blk_ids = output_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1 : num_blks
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
        for blk_elem_index ∈ 1 : num_blk_elems
            element_stress = block_stress[blk_elem_index]
            for point ∈ 1 : num_points
                point_stress = element_stress[point]
                stress_xx[blk_elem_index, point] = point_stress[1]
                stress_yy[blk_elem_index, point] = point_stress[2]
                stress_zz[blk_elem_index, point] = point_stress[3]
                stress_yz[blk_elem_index, point] = point_stress[4]
                stress_xz[blk_elem_index, point] = point_stress[5]
                stress_xy[blk_elem_index, point] = point_stress[6]
            end
        end
        for point ∈ 1 : num_points
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