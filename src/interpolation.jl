function barycentricD2N3(ξ::Vector{Float64})
    N = [1.0 - ξ[1] - ξ[2], ξ[1], ξ[2]]
    dN = [-1 1 0
        -1 0 1] / 1.0
    return N, dN
end

function barycentricD2N3G1()
    w = 0.5 * ones(1, 1)
    N = zeros(3, 1)
    dN = zeros(2, 3, 1)
    ξ = ones(2) / 3.0
    N, dN[:, :, 1] = barycentricD2N3(ξ)
    return N, dN, w
end

function barycentricD2N3G3()
    w = ones(1, 3) / 6.0
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
    w = ones(1, 1) / 6.0
    N = zeros(4, 1)
    dN = zeros(3, 4, 1)
    ξ = 0.25 * ones(3, 1)
    N, dN[:, :, 1] = barycentricD3N4(ξ)
    return N, dN, w
end

function barycentricD3N4G4()
    w = ones(1, 4) / 24.0
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
    w = 2.0 * ones(1, 1)
    N = zeros(2, 1)
    dN = zeros(1, 2, 1)
    N, dN[:, :, 1] = lagrangianD1N2(0.0)
    return N, dN, w
end

function lagrangianD1N2G2()
    w = ones(1, 2)
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
    N = zeros(1, 4)
    dN = zeros(2, 4)
    for i ∈ 1:4
        N[i] = 0.25 * (1.0 + ra[i] * r) * (1.0 + sa[i] * s)
        dN[1, i] = 0.25 * ra[i] * (1 + sa[i] * s)
        dN[2, i] = 0.25 * (1 + ra[i] * r) * sa[i]
    end
    return N, dN
end

function lagrangianD2N4G4()
    w = ones(1, 4)
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
    N = zeros(1, 8)
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
    w = ones(1, 8)
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

#
# Compute isoparametric interpolation functions, their parametric
# derivatives and integration weights.
#
function isoparametric(element_type, num_int)
    msg1 = "Invalid number of integration points: "
    msg2 = " for element type: "
    str_int = string(num_int)
    if element_type == "BAR2"
        if num_int == 1
            N, dN, w = lagrangianD1N2G1()
        elseif num_int == 2
            N, dN, w = lagrangianD1N2G2()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TRI3"
        if num_int == 1
            N, dN, w = barycentricD2N3G1()
        elseif num_int == 3
            N, dN, w = barycentricD2N3G3()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "QUAD4"
        if num_int == 4
            N, dN, w = lagrangianD2N4G4()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "TETRA4"
        if num_int == 1
            N, dN, w = barycentricD3N4G1()
        elseif num_int == 4
            N, dN, w = barycentricD4N4G4()
        else
            error(msg1, num_int, msg2, element_type)
        end
    elseif element_type == "HEX8"
        if num_int == 8
            N, dN, w = lagrangianD3N8G8()
        else
            error(msg1, num_int, msg2, element_type)
        end
    else
        error("Invalid element type: ", element_type)
    end
    return N, dN, w
end

function gradient_operator(dNdX)
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

function get_side_set_nodal_forces(coordinates::Matrix{Float64}, expr::Any, time::Float64)
    _, num_nodes = size(coordinates)
    if num_nodes == 3
        return get_triangle_nodal_forces(coordinates, expr, time)
    elseif num_nodes == 4
        return get_quadrilateral_nodal_forces(coordinates, expr, time)
    else
        error("Unknown side topology with number of nodes: ", num_nodes)
    end
end

function get_triangle_nodal_forces(coordinates3D::Matrix{Float64}, expr::Any, time::Float64)
    A = coordinates3D[:, 1]
    B = coordinates3D[:, 2]
    C = coordinates3D[:, 3]
    centroid = (A + B + C) / 3.0
    coordinates2D = triangle_3D_to_2D(A, B, C)
    Nₚ, dNdξ, elem_weights = barycentricD2N3G1()
    point = 1
    dNdξₚ = dNdξ[:, :, point]
    dXdξ = dNdξₚ * coordinates2D'
    j = det(dXdξ)
    w = elem_weights[point]
    global t = time
    global x = centroid[1]
    global y = centroid[2]
    global z = centroid[3]
    traction = eval(expr)
    nodal_force_component = traction * Nₚ * j * w
    return nodal_force_component
end

function get_quadrilateral_nodal_forces(coordinates3D::Matrix{Float64}, expr::Any, time::Float64)
    A = coordinates3D[:, 1]
    B = coordinates3D[:, 2]
    C = coordinates3D[:, 3]
    D = coordinates3D[:, 4]
    centroid = (A + B + C + D) / 4.0
    coordinates2D = quadrilateral_3D_to_2D(A, B, C, D)
    N, dNdξ, elem_weights = lagrangianD2N4G4()
    g = sqrt(3.0) / 3.0
    global t = time
    nodal_force_component = zeros(4)
    for point ∈ 1:4
        Nₚ = N[:, point]
        dNdξₚ = dNdξ[:, :, point]
        dXdξ = dNdξₚ * coordinates2D'
        j = det(dXdξ)
        w = elem_weights[point]
        point_coordinates3D = g * (coordinates3D[:, point] - centroid)
        global x = point_coordinates3D[1]
        global y = point_coordinates3D[2]
        global z = point_coordinates3D[3]
        traction = eval(expr)
        nodal_force_component += traction * Nₚ * j * w
    end
    return nodal_force_component
end