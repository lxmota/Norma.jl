include("constitutive.jl")
include("interpolation.jl")
include("ics_bcs.jl")

function create_elem_to_blk_map(mesh::PyObject)
    num_blks = mesh.num_blks()
    num_elem = mesh.num_elems()
    elem_to_blk_map = zeros(Int64,num_elem)
    for blk_id ∈ 1:num_blks
        blk_elem_to_blk_map = mesh.get_block_id_map("EX_ELEM_BLOCK", blk_id)
        for elem_id ∈ blk_elem_to_blk_map
            elem_to_blk_map[elem_id] = blk_id
        end
    end
    mesh.elem_to_blk_map = elem_to_blk_map
end

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
    internal_force = zeros(3*num_nodes)
    boundary_traction_force = zeros(3*num_nodes)
    boundary_conditions = Vector{BoundaryCondition}()
    free_dofs = trues(3 * num_nodes)
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
    create_elem_to_blk_map(input_mesh)
    SolidMechanics(input_mesh, materials, reference, current, velocity, acceleration,
        internal_force, boundary_traction_force, boundary_conditions, stress, free_dofs, time, failed)
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
    internal_heat_flux = zeros(num_nodes)
    boundary_heat_flux = zeros(num_nodes)
    boundary_conditions = Vector{BoundaryCondition}()
    free_dofs = trues(num_nodes)
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
    create_elem_to_blk_map(input_mesh)
    HeatConduction(input_mesh, materials, reference, temperature, rate, internal_heat_flux,
    boundary_heat_flux, boundary_conditions, flux, free_dofs, time, failed)
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

function voigt_cauchy_from_stress(_::Solid, P::Matrix{Float64}, F::Matrix{Float64}, J::Float64)
    σ = F * P' ./ J
    return [σ[1, 1], σ[2, 2], σ[3, 3], σ[2, 3], σ[1, 3], σ[1, 2]]
end

function voigt_cauchy_from_stress(_::Linear_Elastic, σ::Matrix{Float64}, _::Matrix{Float64}, _::Float64)
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
    materials = model.materials
    input_mesh = model.mesh
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    energy = 0.0
    internal_force = zeros(num_dof)
    body_force = zeros(num_dof)
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
                F = dxdX
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
    model.internal_force = internal_force
    return energy, internal_force, body_force, stiffness_matrix
end

function evaluate(_::Newmark, model::SolidMechanics)
    materials = model.materials
    input_mesh = model.mesh
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    energy = 0.0
    internal_force = zeros(num_dof)
    body_force = zeros(num_dof)
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
                F = dxdX
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
    model.internal_force = internal_force
    return energy, internal_force, body_force, stiffness_matrix, mass_matrix
end

function get_minimum_edge_length(nodal_coordinates::Matrix{Float64}, edges::Vector{Tuple{Int64,Int64}})
    minimum_edge_length = Inf
    for edge ∈ edges
        node_a = edge[1]
        node_b = edge[2]
        edge_vector = nodal_coordinates[:, node_a] - nodal_coordinates[:, node_b]
        distance = norm(edge_vector)
        minimum_edge_length = min(minimum_edge_length, distance)
    end
    return minimum_edge_length
end

function get_minimum_edge_length(nodal_coordinates::Matrix{Float64}, elem_type::String)
    if elem_type == "TETRA4"
        edges = [(1,2), (1,3), (1,4), (2,3), (3,4), (2,4)]
        return get_minimum_edge_length(nodal_coordinates, edges)
    elseif elem_type == "HEX8"
        edges = [(1,4), (1,5), (4,8), (5,8), (2,3), (2,6), (3,7), (6,7), (1,2), (3,4), (5,6), (7,8)]
        return get_minimum_edge_length(nodal_coordinates, edges)
    else
        error("Invalid element type: ", elem_type)
    end
end

function set_time_step(integrator::CentralDifference, model::SolidMechanics)
    materials = model.materials
    input_mesh = model.mesh
    elem_blk_ids = input_mesh.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    stable_time_step = Inf
    for blk_index ∈ 1:num_blks
        material = materials[blk_index]
        ρ = material.ρ
        M = get_p_wave_modulus(material)
        wave_speed = sqrt(M / ρ)
        minimum_blk_edge_length = Inf
        blk_id = elem_blk_ids[blk_index]
        elem_type = input_mesh.elem_type(blk_id)
        elem_blk_conn, num_blk_elems, num_elem_nodes = input_mesh.get_elem_connectivity(blk_id)
        for blk_elem_index ∈ 1:num_blk_elems
            conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
            node_indices = elem_blk_conn[conn_indices]
            elem_cur_pos = model.current[:, node_indices]
            minimum_elem_edge_length = get_minimum_edge_length(elem_cur_pos, elem_type)
            minimum_blk_edge_length = min(minimum_blk_edge_length, minimum_elem_edge_length)
        end
        blk_stable_time_step = integrator.CFL * minimum_blk_edge_length / wave_speed
        stable_time_step = min(stable_time_step, blk_stable_time_step)
    end
    integrator.stable_time_step = stable_time_step
    if stable_time_step < integrator.user_time_step
        println("Warning: Estimated stable time step: ", stable_time_step, " < provided time step: ", integrator.user_time_step)
    end
    integrator.time_step = min(stable_time_step, integrator.user_time_step)
end

function evaluate(_::CentralDifference, model::SolidMechanics)
    materials = model.materials
    input_mesh = model.mesh
    x, _, _ = input_mesh.get_coords()
    num_nodes = length(x)
    num_dof = 3 * num_nodes
    energy = 0.0
    internal_force = zeros(num_dof)
    body_force = zeros(num_dof)
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
                F = dxdX
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
    model.internal_force = internal_force
    return energy, internal_force, body_force, lumped_mass
end