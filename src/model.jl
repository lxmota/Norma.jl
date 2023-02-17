using SparseArrays

include("constitutive.jl")
include("model_def.jl")
include("time_def.jl")

function SolidMechanics(params::Dict{Any,Any})
    input_mesh = params["input_mesh"]
    model_params = params["model"]
    x, y, z = input_mesh.get_coords()
    num_nodes = length(x)
    reference = Matrix{Float64}(undef, 3, num_nodes)
    current = Matrix{Float64}(undef, 3, num_nodes)
    velocity = Matrix{Float64}(undef, 3, num_nodes)
    acceleration = Matrix{Float64}(undef, 3, num_nodes)
    for node ∈ 1:num_nodes
        reference[:, node] = [x[node], y[node], z[node]]
        current[:, node] = [x[node], y[node], z[node]]
        velocity[:, node] = [0.0, 0.0, 0.0]
        acceleration[:, node] = [0.0, 0.0, 0.0]
    end
    material_params = model_params["material"]
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    if (num_blks_params ≠ num_blks)
        error("number of blocks in mesh ", model_params["mesh"], " (", num_blks,
            ") must be equal to number of blocks in materials file ", model_params["material"],
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
    nodal_dofs = [free::DOF for _ ∈ 1:3*num_nodes]
    stress = Vector{Vector{Vector{Vector{Float64}}}}(undef, num_blks)
    for blk_index ∈ 1:num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        blk_conn = input_mesh.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        block_stress = Vector{Vector{Vector{Float64}}}(undef, num_blk_elems)
        for blk_elem_index ∈ 1:num_blk_elems
            element_stress = Vector{Vector{Float64}}(undef, num_points)
            for point ∈ 1:num_points
                element_stress[point] = zeros(6)
            end
            block_stress[blk_elem_index] = element_stress
        end
        stress[blk_index] = block_stress
    end
    SolidMechanics(params, materials, reference, current, velocity, acceleration, stress, nodal_dofs, time, failed)
end

function HeatConduction(params::Dict{Any,Any})
    input_mesh = params["input_mesh"]
    model_params = params["model"]
    x, y, z = input_mesh.get_coords()
    num_nodes = length(x)
    reference = Matrix{Float64}(undef, 3, num_nodes)
    temperature = Vector{Float64}(undef, num_nodes)
    rate = Vector{Float64}(undef, num_nodes)
    for node ∈ 1:num_nodes
        reference[:, node] = [x[node], y[node], z[node]]
        temperature[node] = 0.0
        rate[node] = 0.0
    end
    material_params = model_params["material"]
    material_blocks = material_params["blocks"]
    num_blks_params = length(material_blocks)
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    if (num_blks_params ≠ num_blks)
        error("number of blocks in mesh ", model_params["mesh"], " (", num_blks,
            ") must be equal to number of blocks in materials file ", model_params["material"],
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
    nodal_dofs = [free::DOF for _ ∈ 1:num_nodes]
    flux = Vector{Vector{Vector{Vector{Float64}}}}(undef, num_blks)
    for blk_index ∈ 1:num_blks
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        blk_conn = input_mesh.get_elem_connectivity(blk_id)
        num_blk_elems = blk_conn[2]
        block_flux = Vector{Vector{Vector{Float64}}}(undef, num_blk_elems)
        for blk_elem_index ∈ 1:num_blk_elems
            element_flux = Vector{Vector{Float64}}(undef, num_points)
            for point ∈ 1:num_points
                element_flux[point] = zeros(3)
            end
            block_flux[blk_elem_index] = element_stress
        end
        flux[blk_index] = block_flux
    end
    HeatConduction(params, materials, reference, temperature, rate, flux, nodal_dofs, time, failed)
end

function create_model(params::Dict{Any,Any})
    model_params = params["model"]
    model_name = model_params["type"]
    if model_name == "solid mechanics"
        return SolidMechanics(params)
    elseif model_name == "heat conduction"
        return HeatConduction(params)
    else
        error("Unknown type of model : ", model_name)
    end
end

function voigt_cauchy_from_stress(_::Any, P::MTTensor, F::MTTensor, J::Float64)
    σ = F * P' ./ J
    return [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]
end

function voigt_cauchy_from_stress(_::Linear_Elastic, σ::MTTensor, _::MTTensor, _::Float64)
    return [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]
end

function assemble(rows::Vector{Int64}, cols::Vector{Int64}, global_stiffness::Vector{Float64}, element_stiffness::Matrix{Float64}, dofs::Vector{Int64})
    num_dofs = length(dofs)
    for i ∈ 1:num_dofs
        I = dofs[i]
        for j ∈ 1:num_dofs
            J = dofs[j]
            push!(rows, I)
            push!(cols, J)
            push!(global_stiffness, element_stiffness[i, j])
        end
    end
end

function assemble(rows::Vector{Int64}, cols::Vector{Int64}, global_stiffness::Vector{Float64}, global_mass::Vector{Float64}, element_stiffness::Matrix{Float64}, element_mass::Matrix{Float64}, dofs::Vector{Int64})
    num_dofs = length(dofs)
    for i ∈ 1:num_dofs
        I = dofs[i]
        for j ∈ 1:num_dofs
            J = dofs[j]
            push!(rows, I)
            push!(cols, J)
            push!(global_mass, element_mass[i, j])
            push!(global_stiffness, element_stiffness[i, j])
        end
    end
end

function evaluate(_::QuasiStatic, model::SolidMechanics)
    params = model.params
    materials = model.materials
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    energy = 0.0
    internal_force = zeros(num_dof)
    external_force = zeros(num_dof)
    rows = Vector{Int64}()
    cols = Vector{Int64}()
    stiffness = Vector{Float64}()
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1:num_blks
        material = materials[blk_index]
        ρ = material.ρ
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        _, dNdξ, elem_weights = isoparametric(elem_type, num_points)
        elem_blk_conn, num_blk_elems, num_elem_nodes = input_mesh.get_elem_connectivity(blk_id)
        num_elem_dofs = 3 * num_elem_nodes
        elem_dofs = zeros(Int64, num_elem_dofs)
        for blk_elem_index ∈ 1:num_blk_elems
            conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
            node_indices = elem_blk_conn[conn_indices]
            elem_ref_pos = model.reference[:, node_indices]
            elem_cur_pos = model.current[:, node_indices]
            element_energy = 0.0
            element_internal_force = zeros(num_elem_dofs)
            element_stiffness = zeros(num_elem_dofs, num_elem_dofs)
            elem_dofs[1:3:num_elem_dofs-2] = 3 .* node_indices .- 2
            elem_dofs[2:3:num_elem_dofs-1] = 3 .* node_indices .- 1
            elem_dofs[3:3:num_elem_dofs] = 3 .* node_indices
            for point ∈ 1:num_points
                dNdξₚ = dNdξ[:, :, point]
                dXdξ = dNdξₚ * elem_ref_pos'
                dxdξ = dNdξₚ * elem_cur_pos'
                dxdX = dXdξ \ dxdξ
                dNdX = dXdξ \ dNdξₚ
                B = gradient_operator(dNdX)
                j = det(dXdξ)
                J = det(dxdX)
                if J ≤ 0.0
                    model.failed = true
                    error("evaluation of model has failed with a non-positive Jacobian")
                    return 0.0, zeros(num_dof), zeros(num_dof), spzeros(num_dof, num_dof)
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
            energy += element_energy
            internal_force[elem_dofs] += element_internal_force
            assemble(rows, cols, stiffness, element_stiffness, elem_dofs)
        end
    end
    stiffness_matrix = sparse(rows, cols, stiffness)
    return energy, internal_force, external_force, stiffness_matrix
end

function evaluate(_::Newmark, model::SolidMechanics)
    params = model.params
    materials = model.materials
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    energy = 0.0
    internal_force = zeros(num_dof)
    external_force = zeros(num_dof)
    rows = Vector{Int64}()
    cols = Vector{Int64}()
    stiffness = Vector{Float64}()
    mass = Vector{Float64}()
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1:num_blks
        material = materials[blk_index]
        ρ = material.ρ
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        N, dNdξ, elem_weights = isoparametric(elem_type, num_points)
        elem_blk_conn, num_blk_elems, num_elem_nodes = input_mesh.get_elem_connectivity(blk_id)
        num_elem_dofs = 3 * num_elem_nodes
        elem_dofs = zeros(Int64, num_elem_dofs)
        for blk_elem_index ∈ 1:num_blk_elems
            conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
            node_indices = elem_blk_conn[conn_indices]
            elem_ref_pos = model.reference[:, node_indices]
            elem_cur_pos = model.current[:, node_indices]
            element_energy = 0.0
            element_internal_force = zeros(num_elem_dofs)
            element_stiffness = zeros(num_elem_dofs, num_elem_dofs)
            element_mass = zeros(num_elem_dofs, num_elem_dofs)
            index_x = 1:3:num_elem_dofs.-2
            index_y = index_x .+ 1
            index_z = index_x .+ 2
            elem_dofs[index_x] = 3 .* node_indices .- 2
            elem_dofs[index_y] = 3 .* node_indices .- 1
            elem_dofs[index_z] = 3 .* node_indices
            for point ∈ 1:num_points
                dNdξₚ = dNdξ[:, :, point]
                dXdξ = dNdξₚ * elem_ref_pos'
                dxdξ = dNdξₚ * elem_cur_pos'
                dxdX = dXdξ \ dxdξ
                dNdX = dXdξ \ dNdξₚ
                B = gradient_operator(dNdX)
                j = det(dXdξ)
                J = det(dxdX)
                if J ≤ 0.0
                    model.failed = true
                    error("evaluation of model has failed with a non-positive Jacobian")
                    return 0.0, zeros(num_dof), zeros(num_dof), spzeros(num_dof, num_dof), spzeros(num_dof, num_dof)
                end
                F = MTTensor(dxdX)
                W, P, A = constitutive(material, F)
                stress = P[1:9]
                moduli = second_from_fourth(A)
                w = elem_weights[point]
                element_energy += W * j * w
                element_internal_force += B' * stress * j * w
                element_stiffness += B' * moduli * B * j * w
                reduced_mass = N[:, point] * N[:, point]' * ρ * j * w
                element_mass[index_x, index_x] += reduced_mass
                element_mass[index_y, index_y] += reduced_mass
                element_mass[index_z, index_z] += reduced_mass
                voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
                model.stress[blk_index][blk_elem_index][point] = voigt_cauchy
            end
            energy += element_energy
            internal_force[elem_dofs] += element_internal_force
            assemble(rows, cols, stiffness, mass, element_stiffness, element_mass, elem_dofs)
        end
    end
    stiffness_matrix = sparse(rows, cols, stiffness)
    mass_matrix = sparse(rows, cols, mass)
    return energy, internal_force, external_force, stiffness_matrix, mass_matrix
end

function get_minimum_edge_length(nodal_coordinates::Matrix{Float64})
    _, num_nodes = size(nodal_coordinates)
    minimum_edge_length = Inf
    origin = nodal_coordinates[:, 1]
    for node = 2:num_nodes
        edge = nodal_coordinates[:, node] - origin
        distance = sqrt(edge' * edge)
        minimum_edge_length = min(minimum_edge_length, distance)
    end
    return minimum_edge_length
end

function set_time_step(integrator::CentralDifference, model::SolidMechanics)
    params = model.params
    materials = model.materials
    input_mesh = params["input_mesh"]
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    stable_time_step = Inf
    for blk_index ∈ 1:num_blks
        material = materials[blk_index]
        ρ = material.ρ
        M = estimate_p_wave_modulus(material)
        wave_speed = sqrt(M / ρ)
        minimum_blk_edge_length = Inf
        blk_id = elem_blk_ids[blk_index]
        elem_blk_conn, num_blk_elems, num_elem_nodes = input_mesh.get_elem_connectivity(blk_id)
        for blk_elem_index ∈ 1:num_blk_elems
            conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
            node_indices = elem_blk_conn[conn_indices]
            elem_cur_pos = model.current[:, node_indices]
            minimum_elem_edge_length = get_minimum_edge_length(elem_cur_pos)
            minimum_blk_edge_length = min(minimum_blk_edge_length, minimum_elem_edge_length)
        end
        blk_stable_time_step = minimum_blk_edge_length / wave_speed
        stable_time_step = integrator.CFL * min(stable_time_step, blk_stable_time_step)
    end
    integrator.stable_time_step = stable_time_step
    if stable_time_step < integrator.user_time_step
        println("Warning: Estimated stable time step: ", stable_time_step, " < provided time step: ", integrator.user_time_step)
    end
    integrator.time_step = min(stable_time_step, integrator.user_time_step)
end

function evaluate(_::CentralDifference, model::SolidMechanics)
    params = model.params
    materials = model.materials
    input_mesh = params["input_mesh"]
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    energy = 0.0
    internal_force = zeros(num_dof)
    external_force = zeros(num_dof)
    lumped_mass = zeros(num_dof)
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for blk_index ∈ 1:num_blks
        material = materials[blk_index]
        ρ = material.ρ
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        num_points = default_num_int_pts(elem_type)
        N, dNdξ, elem_weights = isoparametric(elem_type, num_points)
        elem_blk_conn, num_blk_elems, num_elem_nodes = input_mesh.get_elem_connectivity(blk_id)
        num_elem_dofs = 3 * num_elem_nodes
        elem_dofs = zeros(Int64, num_elem_dofs)
        for blk_elem_index ∈ 1:num_blk_elems
            conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
            node_indices = elem_blk_conn[conn_indices]
            elem_ref_pos = model.reference[:, node_indices]
            elem_cur_pos = model.current[:, node_indices]
            element_energy = 0.0
            element_internal_force = zeros(num_elem_dofs)
            element_lumped_mass = zeros(num_elem_dofs)
            index_x = 1:3:num_elem_dofs.-2
            index_y = index_x .+ 1
            index_z = index_x .+ 2
            elem_dofs[index_x] = 3 .* node_indices .- 2
            elem_dofs[index_y] = 3 .* node_indices .- 1
            elem_dofs[index_z] = 3 .* node_indices
            for point ∈ 1:num_points
                dNdξₚ = dNdξ[:, :, point]
                dXdξ = dNdξₚ * elem_ref_pos'
                dxdξ = dNdξₚ * elem_cur_pos'
                dxdX = dXdξ \ dxdξ
                dNdX = dXdξ \ dNdξₚ
                B = gradient_operator(dNdX)
                j = det(dXdξ)
                J = det(dxdX)
                if J ≤ 0.0
                    model.failed = true
                    error("evaluation of model has failed with a non-positive Jacobian")
                    return 0.0, zeros(num_dof), zeros(num_dof), zeros(num_dof)
                end
                F = MTTensor(dxdX)
                W, P, _ = constitutive(material, F)
                stress = P[1:9]
                w = elem_weights[point]
                element_energy += W * j * w
                element_internal_force += B' * stress * j * w
                reduced_mass = N[:, point] * N[:, point]' * ρ * j * w
                reduced_lumped_mass = sum(reduced_mass, dims=2)
                element_lumped_mass[index_x] += reduced_lumped_mass
                element_lumped_mass[index_y] += reduced_lumped_mass
                element_lumped_mass[index_z] += reduced_lumped_mass
                voigt_cauchy = voigt_cauchy_from_stress(material, P, F, J)
                model.stress[blk_index][blk_elem_index][point] = voigt_cauchy
            end
            energy += element_energy
            internal_force[elem_dofs] += element_internal_force
            lumped_mass[elem_dofs] += element_lumped_mass
        end
    end
    return energy, internal_force, external_force, lumped_mass
end

function node_set_id_from_name(node_set_name::String, mesh::PyObject)
    node_set_names = mesh.get_node_set_names()
    num_names = length(node_set_names)
    node_set_index = 0
    for index ∈ 1:num_names
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
    if haskey(params, "boundary conditions") == false
        return
    end
    reference = model.reference
    current = model.current
    input_mesh = params["input_mesh"]
    global t = model.time
    _, num_nodes = size(reference)
    nodal_dofs = [free::DOF for _ ∈ 1:3*num_nodes]
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
                # function_str is an arbitrary function of t, x, y, z in the input file
                bc_expr = Meta.parse(function_str)
                for node_index ∈ node_set_node_indices
                    global x = reference[1, node_index]
                    global y = reference[2, node_index]
                    global z = reference[3, node_index]
                    bc_val = eval(bc_expr)
                    dof_index = 3 * (node_index - 1) + offset
                    current[offset, node_index] = reference[offset, node_index] + bc_val
                    nodal_dofs[dof_index] = Dirichlet::DOF
                end
            elseif bc_type == "Schwarz"
            elseif bc_type == "Neumann"
            end
        end
    end
    model.nodal_dofs = nodal_dofs
end

function apply_ics(model::SolidMechanics)
    params = model.params
    if haskey(params, "initial conditions") == false
        return
    end
    reference = model.reference
    current = model.current
    velocity = model.velocity
    input_mesh = params["input_mesh"]
    global t = model.time
    ic_params = params["initial conditions"]
    for (ic_type, ic_type_params) ∈ ic_params
        for ic ∈ ic_type_params
            node_set_name = ic["node set"]
            function_str = ic["function"]
            component = ic["component"]
            offset = component_offset_from_string(component)
            node_set_id = node_set_id_from_name(node_set_name, input_mesh)
            node_set_node_indices = input_mesh.get_node_set_nodes(node_set_id)
            # function_str is an arbitrary function of t, x, y, z in the input file
            ic_expr = Meta.parse(function_str)
            for node_index ∈ node_set_node_indices
                global x = reference[1, node_index]
                global y = reference[2, node_index]
                global z = reference[3, node_index]
                ic_val = eval(ic_expr)
                if ic_type == "displacement"
                    current[offset, node_index] = reference[offset, node_index] + ic_val
                elseif ic_type == "velocity"
                    velocity[offset, node_index] = ic_val
                end
            end
        end
    end
end

function apply_bcs(model::HeatConduction)
    params = model.params
    if haskey(params, "boundary conditions") == false
        return
    end
    temperature = model.temperature
    input_mesh = params["input_mesh"]
    global t = model.time
    xc, yc, zc = input_mesh.get_coords()
    num_nodes = length(xc)
    nodal_dofs = [free::DOF for _ ∈ 1:num_nodes]
    bc_params = params["boundary conditions"]
    for (bc_type, bc_type_params) ∈ bc_params
        for bc ∈ bc_type_params
            if bc_type == "Dirichlet"
                node_set_name = bc["node set"]
                function_str = bc["function"]
                node_set_id = node_set_id_from_name(node_set_name, input_mesh)
                node_set_node_indices = input_mesh.get_node_set_nodes(node_set_id)
                # function_str is an arbitrary function of t, x, y, z in the input file
                bc_expr = Meta.parse(function_str)
                for node_index ∈ node_set_node_indices
                    global x = xc[node_index]
                    global y = yc[node_index]
                    global z = zc[node_index]
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