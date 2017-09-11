#
#
function barycentricD2N3(ξ)
    Na  = [1 - ξ[1] - ξ[2],
           ξ[1],
           ξ[2]]
    dNa = [-1  1  0;
           -1  0  1]
    return Na, dNa
end

#
#
function barycentricD2N3G1()
    w = 0.5 * ones(1, 1)
    Na = zeros(3, 1)
    dNa = zeros(2, 3, 1)
    ξ = ones(2, 1) / 3.0
    Na, dNa[:, :, 1] = barycentricD2N3(ξ)
    return Na, dNa, w
end

#
#
function barycentricD2N3G3()
    w = ones(1, 3) / 6.0
    Na = zeros(3, 3)
    dNa = zeros(2, 3, 3)
    ξ = [4  1  1  1;
         1  4  1  1;
         1  1  4  1] / 6.0
    for i = 1 : 3
        Na[:, i], dNa[:, :, i] = barycentricD2N3(ξ[:, i])
    end
    return Na, dNa, w
end

#
#
function barycentricD3N4(ξ)
    Na  = [1 - ξ[1] - ξ[2] - ξ[3],
           ξ[1],
           ξ[2],
           ξ[3]]
    dNa = [-1  1  0  0;
           -1  0  1  0;
           -1  0  0  1]
    return Na, dNa
end

#
#
function barycentricD3N4G1()
    w = ones(1, 1) / 6.0
    Na = zeros(4, 1)
    dNa = zeros(3, 4, 1)
    ξ = 0.25 * ones(3, 1)
    Na, dNa[:, :, 1] = barycentricD3N4(ξ)
    return Na, dNa, w
end

#
#
function barycentricD3N4G4()
    w = ones(1, 4) / 24.0
    Na = zeros(4, 4)
    dNa = zeros(3, 4, 4)
    s = sqrt(5.0)
    a = 5.0 + 3.0 * s
    b = 5.0 - s
    ξ = [a  b  b  b;
         b  a  b  b;
         b  b  a  b] / 20.0
    for i = 1 : 4
        Na[:, i], dNa[:, :, i] = barycentricD3N4(ξ[:, i])
    end
    return Na, dNa, w
end

#
#
function lagrangianD1N2(ξ)
    Na  = [0.5 * (1.0 - ξ), 0.5 * (1.0 + ξ)]  
    dNa = [-0.5, 0.5]
    return Na, dNa
end

#
#
function lagrangianD1N2G1()
    w = 2.0 * ones(1, 1)
    Na = zeros(2, 1)
    dNa = zeros(1, 2, 1)
    Na, dNa[:, :, 1] = lagrangianD1N2(0.0)
    return Na, dNa, w
end

#
#
function lagrangianD1N2G2()
    w = ones(1, 2)
    Na = zeros(2, 2)
    dNa = zeros(1, 2, 2)
    g = sqrt(3.0) / 3.0
    ξ = [-g, g]
    for i = 1 : 2
        Na[:, i], dNa[:, :, i] = lagrangianD1N2(ξ[i])
    end
    return Na, dNa, w
end

#
#
function lagrangianD2N4(ξ)
    r = ξ[1]
    s = ξ[2]
    ra = [-1  1  1 -1]
    sa = [-1 -1  1  1]
    Na  = zeros(1, 4)
    dNa = zeros(2, 4)
    for i = 1 : 4
        Na[i] = 0.25 * (1 + ra[i] * r) * (1 + sa[i] * s)
        dNa[1, i] = 0.25 * ra[i] * (1 + sa[i] * s)
        dNa[2, i] = 0.25 * (1 + ra[i] * r) * sa[i]
    end
    return Na, dNa
end
function lagrangianD2N4G4()
    w = ones(1, 4)
    Na = zeros(4, 4)
    dNa = zeros(2, 4, 4)  
    g = sqrt(3.0) / 3.0
    ξ = g * [-1  1  1 -1;
             -1 -1  1  1]
    for i = 1 : 4
        Na[:, i], dNa[:, :, i] = lagrangianD2N4(ξ[:, i])
    end
    return Na, dNa, w
end

#
#
function lagrangianD3N8(ξ)
    r = ξ[1]
    s = ξ[2]
    t = ξ[3]
    ra = [-1  1  1 -1 -1  1  1 -1]
    sa = [-1 -1  1  1 -1 -1  1  1]
    ta = [-1 -1 -1 -1  1  1  1  1]
    Na  = zeros(1, 8)
    dNa = zeros(3, 8)
    for i = 1 : 8
        Na[i] = 0.125 * (1 + ra[i] * r) * (1 + sa[i] * s) * (1 + ta[i] * t)
        dNa[1, i] = 0.125 * ra[i] * (1 + sa[i] * s) * (1 + ta[i] * t)
        dNa[2, i] = 0.125 * (1 + ra[i] * r) * sa[i] * (1 + ta[i] * t)
        dNa[3, i] = 0.125 * (1 + ra[i] * r) * (1 + sa[i] * s) * ta[i]
    end
    return Na, dNa
end

#
#
function lagrangianD3N8G8()
    w = ones(1, 8)
    Na = zeros(8, 8)
    dNa = zeros(3, 8, 8)
    g = sqrt(3.0) / 3.0
    ξ = g * [-1  1  1 -1 -1  1  1 -1;
             -1 -1  1  1 -1 -1  1  1;
             -1 -1 -1 -1  1  1  1  1]
    for i = 1 : 8
        Na[:, i], dNa[:, :, i] = lagrangianD3N8(ξ[:, i])
    end
    return Na, dNa, w
end

#
# Compute isoparametric interpolation functions, their parametric
# derivatives and integration weights.
#
function isoparametric(element_type, num_int)
    msg = "Invalid number of integration points for element: "
    str_int = string(num_int)
    if element_type == "BAR"
        if num_int == 1
            Na, dNa, w = lagrangianD1N2G1()
        elseif num_int == 2
            Na, dNa, w = lagrangianD1N2G2()
        else
            error("$msg $str_int $element_type")
        end
    elseif element_type == "TRI3"
        if num_int == 1
            Na, dNa, w = barycentricD2N3G1()
        elseif num_int == 3
            Na, dNa, w = barycentricD2N3G3()
        else
            error("$msg $str_int $element_type")
        end
    elseif element_type == "QUAD4"
        if num_int == 4
            Na, dNa, w = lagrangianD2N4G4()
        else
            error("$msg $str_int $element_type")
        end
    elseif element_type == "TETRA4"
        if num_int == 1
            Na, dNa, w = barycentricD3N4G1()
        elseif num_int == 4
            Na, dNa, w = barycentricD4N4G4()
        else
            error("$msg $str_int $element_type")
        end
    elseif element_type == "HEX8"
        if num_int == 8
            Na, dNa, w = lagrangianD3N8G8()
        else
            error("$msg $str_int $element_type")
        end
    else
        error("invalid element type: $element_type")
    end
    return Na, dNa, w
end
