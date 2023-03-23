function barycentricD2N3(ξ::Vector{Float64})
    N = [1.0 - ξ[1] - ξ[2], ξ[1], ξ[2]]
    dN = [-1 1 0
        -1 0 1] / 1.0
    return N, dN
end

function barycentricD2N3G1()
    w = 0.5 * ones(1)
    N = zeros(3, 1)
    dN = zeros(2, 3, 1)
    ξ = ones(2) / 3.0
    N, dN[:, :, 1] = barycentricD2N3(ξ)
    return N, dN, w
end

function barycentricD2N3G3()
    w = ones(3) / 6.0
    N = zeros(3, 3)
    dN = zeros(2, 3, 3)
    ξ = [4 1 1 1
        1 4 1 1
        1 1 4 1] / 6.0
    for i ∈ 1:3
        N[:, i], dN[:, :, i] = barycentricD2N3(ξ[:, i])
    end
    return N, dN, w
end

function barycentricD3N4(ξ::Vector{Float64})
    N = [1.0 - ξ[1] - ξ[2] - ξ[3],
        ξ[1],
        ξ[2],
        ξ[3]]
    dN = [-1 1 0 0
        -1 0 1 0
        -1 0 0 1] / 1.0
    return N, dN
end

function barycentricD3N4G1()
    w = ones(1) / 6.0
    N = zeros(4, 1)
    dN = zeros(3, 4, 1)
    ξ = 0.25 * ones(3, 1)
    N, dN[:, :, 1] = barycentricD3N4(ξ)
    return N, dN, w
end

function barycentricD3N4G4()
    w = ones(4) / 24.0
    N = zeros(4, 4)
    dN = zeros(3, 4, 4)
    s = sqrt(5.0)
    a = 5.0 + 3.0 * s
    b = 5.0 - s
    ξ = [a b b b
        b a b b
        b b a b] / 20.0
    for i ∈ 1:4
        N[:, i], dN[:, :, i] = barycentricD3N4(ξ[:, i])
    end
    return N, dN, w
end

function lagrangianD1N2(ξ::Float64)
    N = [0.5 * (1.0 - ξ), 0.5 * (1.0 + ξ)]
    dN = [-0.5, 0.5]
    return N, dN
end

function lagrangianD1N2G1()
    w = 2.0 * ones(1)
    N = zeros(2, 1)
    dN = zeros(1, 2, 1)
    N, dN[:, :, 1] = lagrangianD1N2(0.0)
    return N, dN, w
end

function lagrangianD1N2G2()
    w = ones(2)
    N = zeros(2, 2)
    dN = zeros(1, 2, 2)
    g = sqrt(3.0) / 3.0
    ξ = [-g, g]
    for i ∈ 1:2
        N[:, i], dN[:, :, i] = lagrangianD1N2(ξ[i])
    end
    return N, dN, w
end

function lagrangianD2N4(ξ::Vector{Float64})
    r = ξ[1]
    s = ξ[2]
    ra = [-1 1 1 -1] / 1.0
    sa = [-1 -1 1 1] / 1.0
    N = zeros(4)
    dN = zeros(2, 4)
    for i ∈ 1:4
        N[i] = 0.25 * (1.0 + ra[i] * r) * (1.0 + sa[i] * s)
        dN[1, i] = 0.25 * ra[i] * (1 + sa[i] * s)
        dN[2, i] = 0.25 * (1 + ra[i] * r) * sa[i]
    end
    return N, dN
end

function lagrangianD2N4G4()
    w = ones(4)
    N = zeros(4, 4)
    dN = zeros(2, 4, 4)
    g = sqrt(3.0) / 3.0
    ξ = g * [-1 1 1 -1
        -1 -1 1 1]
    for i ∈ 1:4
        N[:, i], dN[:, :, i] = lagrangianD2N4(ξ[:, i])
    end
    return N, dN, w
end

function lagrangianD3N8(ξ::Vector{Float64})
    r = ξ[1]
    s = ξ[2]
    t = ξ[3]
    ra = [-1 1 1 -1 -1 1 1 -1] / 1.0
    sa = [-1 -1 1 1 -1 -1 1 1] / 1.0
    ta = [-1 -1 -1 -1 1 1 1 1] / 1.0
    N = zeros(8)
    dN = zeros(3, 8)
    for i ∈ 1:8
        N[i] = 0.125 * (1.0 + ra[i] * r) * (1.0 + sa[i] * s) * (1.0 + ta[i] * t)
        dN[1, i] = 0.125 * ra[i] * (1.0 + sa[i] * s) * (1.0 + ta[i] * t)
        dN[2, i] = 0.125 * (1.0 + ra[i] * r) * sa[i] * (1.0 + ta[i] * t)
        dN[3, i] = 0.125 * (1.0 + ra[i] * r) * (1.0 + sa[i] * s) * ta[i]
    end
    return N, dN
end

function lagrangianD3N8G8()
    w = ones(8)
    N = zeros(8, 8)
    dN = zeros(3, 8, 8)
    g = sqrt(3.0) / 3.0
    ξ = g * [-1 1 1 -1 -1 1 1 -1
        -1 -1 1 1 -1 -1 1 1
        -1 -1 -1 -1 1 1 1 1]
    for i ∈ 1:8
        N[:, i], dN[:, :, i] = lagrangianD3N8(ξ[:, i])
    end
    return N, dN, w
end

function default_num_int_pts(element_type)
    if element_type == "BAR2"
        return 1
    elseif element_type == "TRI3"
        return 1
    elseif element_type == "QUAD4"
        return 4
    elseif element_type == "TETRA4"
        return 1
    elseif element_type == "HEX8"
        return 8
    else
        error("Invalid element type: ", element_type)
    end
end

function get_element_type(dim::Int64, num_nodes::Int64)
    if dim == 1 && num_nodes == 2
        return "BAR2"
    elseif dim == 2 && num_nodes == 3
        return "TRI3"
    elseif dim == 2 && num_nodes == 4
        return "QUAD4"
    elseif dim == 3 && num_nodes == 4
        return "TETRA4"
    elseif dim == 3 && num_nodes == 8
        return "HEX8"
    else
        error("Invalid dimension : ", dim, " and number of nodes : ", num_nodes)
    end
end

#
# Compute isoparametric interpolation functions, their parametric
# derivatives and integration weights.
#
function isoparametric(element_type::String, num_int::Int64)
    msg1 = "Invalid number of integration points: "
    msg2 = " for element type: "
    if element_type == "BAR2"
        if num_int == 1
            return lagrangianD1N2G1()
        elseif num_int == 2
            return lagrangianD1N2G2()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TRI3"
        if num_int == 1
            return barycentricD2N3G1()
        elseif num_int == 3
            return barycentricD2N3G3()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "QUAD4"
        if num_int == 4
            return lagrangianD2N4G4()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TETRA4"
        if num_int == 1
            return barycentricD3N4G1()
        elseif num_int == 4
            return barycentricD4N4G4()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "HEX8"
        if num_int == 8
            return lagrangianD3N8G8()
        else
            error(msg1, num_int, msg2, element_type)
        end
    else
        error("Invalid element type: ", element_type)
    end
end

function gradient_operator(dNdX::Matrix{Float64})
    dim, nen = size(dNdX)
    B = zeros(dim * dim, nen * dim)
    for i ∈ 1:dim
        for j ∈ 1:dim
            p = dim * (i - 1) + j
            for a ∈ 1:nen
                for k ∈ 1:dim
                    q = dim * (a - 1) + k
                    B[p, q] = I[i, k] * dNdX[j, a]
                end
            end
        end
    end
    return B
end

function surface_3D_to_2D(vertices::Matrix{Float64})
    _, num_nodes = size(vertices)
    if num_nodes == 3
        return triangle_3D_to_2D(vertices[:, 1], vertices[:, 2], vertices[:, 3])
    elseif num_nodes == 4
        return quadrilateral_3D_to_2D(vertices[:, 1], vertices[:, 2], vertices[:, 3], vertices[:, 4])
    else
        error("Unknown surface type with number of nodes : ", num_nodes)
    end
end

function triangle_3D_to_2D(A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64})
    AB = B - A
    Bx = norm(AB)
    AC = C - A
    AC2 = AC ⋅ AC
    n = AB / Bx
    Cx = AC ⋅ n
    Cy = sqrt(AC2 - Cx * Cx)
    coordinates = zeros(2, 3)
    coordinates[1, 2] = Bx
    coordinates[1, 3] = Cx
    coordinates[2, 3] = Cy
    return coordinates
end

function quadrilateral_3D_to_2D(A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64}, D::Vector{Float64})
    AB = A - B
    AC = A - C
    AD = A - D
    V = cross(AB, AC)
    v = V / norm(V)
    u = AD ⋅ v
    if abs(u) / norm(AD) >= 1.0e-04
        error("Curved quadrilaterals are not supported")
    end
    Bx = norm(AB)
    AC2 = AC ⋅ AC
    AD2 = AD ⋅ AD
    n = AB / Bx
    Cx = AC ⋅ n
    Cy = sqrt(AC2 - Cx * Cx)
    Dx = AD ⋅ n
    Dy = sqrt(AD2 - Dx * Dx)
    coordinates = zeros(2, 4)
    coordinates[1, 2] = Bx
    coordinates[1, 3] = Cx
    coordinates[1, 4] = Dx
    coordinates[2, 3] = Cy
    coordinates[2, 4] = Dy
    return coordinates
end

using Symbolics
@variables t, x, y, z

function get_side_set_nodal_forces(nodal_coord::Matrix{Float64}, traction_num::Num, time::Float64)
    _, num_side_nodes = size(nodal_coord)
    if num_side_nodes == 3
        A = nodal_coord[:, 1]
        B = nodal_coord[:, 2]
        C = nodal_coord[:, 3]
        centroid = (A + B + C) / 3.0
        two_dim_coord = triangle_3D_to_2D(A, B, C)
        g = 0.0
    elseif num_side_nodes == 4
        A = nodal_coord[:, 1]
        B = nodal_coord[:, 2]
        C = nodal_coord[:, 3]
        D = nodal_coord[:, 4]
        centroid = (A + B + C + D) / 4.0
        two_dim_coord = quadrilateral_3D_to_2D(A, B, C, D)
        g = sqrt(3.0) / 3.0
    else
        error("Unknown side topology with number of nodes: ", num_side_nodes)
    end
    element_type = get_element_type(2, num_side_nodes)
    num_int_points = default_num_int_pts(element_type)
    N, dNdξ, w = isoparametric(element_type, num_int_points)
    nodal_force_component = zeros(num_side_nodes)
    for point ∈ 1:num_int_points
        Nₚ = N[:, point]
        dNdξₚ = dNdξ[:, :, point]
        dXdξ = dNdξₚ * two_dim_coord'
        j = det(dXdξ)
        wₚ = w[point]
        point_coord = g * nodal_coord[:, point] + (1.0 - g) * centroid
        values = Dict(t=>time, x=>point_coord[1], y=>point_coord[2], z=>point_coord[3])
        traction_sym = substitute(traction_num, values)
        traction_val = extract_value(traction_sym)
        nodal_force_component += traction_val * Nₚ * j * wₚ
    end
    return nodal_force_component
end

function map_to_parametric(element_type::String, vertices::Matrix{Float64}, point::Vector{Float64})
    tol = 1.0e-08
    dim = length(point)
    ξ = zeros(dim)
    hessian = zeros(dim, dim)
    while true
        N, dN = interpolate(element_type, ξ)
        trial_point = vertices * N
        residual = trial_point - point
        hessian = vertices * dN'
        δ = - hessian \ residual
        ξ = ξ + δ
        error = norm(δ)
        if error <= tol
          break
        end
    end
    return ξ
end

function interpolate(element_type::String, ξ::Vector{Float64})
    if element_type == "BAR2"
        return lagrangianD1N2(ξ)
    elseif element_type == "TRI3"
        return barycentricD2N3(ξ)
    elseif element_type == "QUAD4"
        return lagrangianD2N4(ξ)
    elseif element_type == "TETRA4"
        return barycentricD3N4(ξ)
    elseif element_type == "HEX8"
        return lagrangianD3N8(ξ)
    else
        error("Invalid element type: ", element_type)
    end
end

function is_inside_parametric(element_type::String, ξ::Vector{Float64})
    tol = 1.0e-06
    # Shrink slightly so that if ξ is nearly inside it still counts as inside
    ξ *= (1.0 - tol)
    if element_type == "BAR2"
        return -1.0 ≤ ξ ≤ 1.0
    elseif element_type == "TRI3"
        return reduce(*, zeros(2) .≤ ξ .≤ ones(2))
    elseif element_type == "QUAD4"
        return reduce(*, -ones(2) .≤ ξ .≤ ones(2))
    elseif element_type == "TETRA4"
        return reduce(*, zeros(3) .≤ ξ .≤ ones(3))
    elseif element_type == "HEX8"
        return reduce(*, -ones(3) .≤ ξ .≤ ones(3))
    else
        error("Invalid element type: ", element_type)
    end
end

function is_inside(element_type::String, vertices::Matrix{Float64}, point::Vector{Float64})
    ξ = map_to_parametric(element_type, vertices, point)
    return is_inside_parametric(element_type, ξ)
end

function find_and_project(point::Vector{Float64}, coupled_mesh::PyObject, coupled_block_id::Int64, coupled_side_set_id::Int64, model::SolidMechanics)
    elem_blk_conn, num_blk_elems, num_elem_nodes = coupled_mesh.get_elem_connectivity(coupled_block_id)
    elem_type = coupled_mesh.elem_type(blk_id)
    point_new = point
    for blk_elem_index ∈ 1:num_blk_elems
        conn_indices = (blk_elem_index-1)*num_elem_nodes+1:blk_elem_index*num_elem_nodes
        node_indices = elem_blk_conn[conn_indices]
        vertices = model.current[:, node_indices]
        inside = is_inside(elem_type, vertices, point)
        #if a point is inside an element, we will move it on the contact side
        if inside == true
            #call a function wich projects the point onto the contact boundary
            point_new = project_onto_contact_surface(point, coupled_side_set_id, coupled_mesh, model)
            break
        end        
    end
    if inside == false
        error("Point : ", point, " not in contact")
    end
    return point_new, node_indices
end

function project_onto_contact_surface(point::Vector{Float64}, coupled_side_set_id::Int64, coupled_mesh::PyObject, model::SolidMechanics)
    #we assume that we know the contact surfaces in advance 
    num_nodes_per_side, side_set_node_indices = coupled_mesh.get_side_set_node_list(coupled_side_set_id)
    coupled_ss_node_index = 1
    minimum_distance = Inf
    point_new = point
    for coupled_side ∈ num_nodes_per_side
        coupled_side_nodes = side_set_node_indices[coupled_ss_node_index:coupled_ss_node_index+coupled_side-1]
        coupled_side_coordinates = model.current[:, coupled_side_nodes]
        #plane 
        coordinates_A = coupled_side_coordinates[:, 1]
        coordinates_B = coupled_side_coordinates[:, 2]
        coordinates_C = coupled_side_coordinates[:, 3]
        BA = coordinates_B - coordinates_A
        CA = coordinates_C - coordinates_A
        N = cross(BA, CA)
        n = N / norm(N)
        distance = (point - coordinates_A) ⋅ n
        #store the new point if the distance is min
        if abs(distance) < minimum_distance
            point_new = point - distance * n   
        end    
        minimum_distance = min(minimum_distance, abs(distance))
        coupled_ss_node_index += coupled_side
    end
    #point_new: new points position
    return point_new
end

function get_projection_square_matrix(mesh::PyObject, side_set_id::Int64)
    num_nodes_per_side, side_set_node_indices = coupled_mesh.get_side_set_node_list(coupled_side_set_id)
    num_nodes = length(side_set_node_indices)
    projection_square_matrix = zeros(num_nodes, num_nodes)
    coupled_ss_node_index = 1
    for side ∈ num_nodes_per_side
        side_nodes = side_set_node_indices[coupled_ss_node_index:coupled_ss_node_index+side-1]
        side_coordinates = model.reference[:, side_nodes]
        two_dim_coord = surface_3D_to_2D(side_coordinates)
        num_side_nodes = length(side)
        element_type = get_element_type(2, num_side_nodes)
        num_int_points = default_num_int_pts(element_type)
        N, dNdξ, w = isoparametric(element_type, num_int_points)
        side_matrix = zeros(num_side_nodes, num_side_nodes)
        for point ∈ 1:num_int_points
            Nₚ = N[:, point]
            dNdξₚ = dNdξ[:, :, point]
            dXdξ = dNdξₚ * two_dim_coord'
            j = det(dXdξ)
            wₚ = w[point]
            side_matrix += Nₚ * Nₚ' * j * wₚ
        end
        coupled_ss_node_index += coupled_side
    end
end

function interpolate(tᵃ::Float64, tᵇ::Float64, xᵃ::Vector{Float64}, xᵇ::Vector{Float64}, t::Float64)
    Δt = tᵇ - tᵃ
    if Δt == 0.0
        return 0.5 * (xᵃ + xᵇ)
    end
    p = (tᵇ - t) / Δt
    q = 1.0 - p
    return p * xᵃ + q * xᵇ
end