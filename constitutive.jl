#
# κ λ μ ν ρ
#
struct SaintVenant_Kirchhoff
    E
    ν
    λ
    μ
    ρ
    function SaintVenant_Kirchhoff(E, ν, ρ)
        λ = E * ν / (1.0 + ν) / (1.0 - 2.0 * ν) 
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, λ, μ, ρ)
    end
end

#
#
struct Neohookean
    E
    ν
    κ
    μ
    ρ
    function Neohookean(E, ν, ρ)
        κ = E / (1.0 - 2.0 * ν) / 3.0
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, κ, μ, ρ)
    end
end

#
#
function constitutive(mat::SaintVenant_Kirchhoff, C)
    dim = size(C, 1)
    I = eye(dim)
    E = 0.5 * (C - I)
    λ = mat.λ
    μ = mat.μ
    trE = trace(E)
    W = 0.5 * λ * trE * trE + μ * trace(E * E)
    S = λ * trE * I + 2.0 * μ * E
    CC = zeros(dim, dim, dim, dim)
    for i = 1 : dim
        for j = 1 : dim
            for k = 1 : dim
                for l = 1 : dim
                    DijDkl = I[i, j] * I[k, l]
                    DikDjl = I[i, k] * I[j, l]
                    DilDjk = I[i, l] * I[j, k]
                    CC[i, j, k, l] = λ * DijDkl + μ * (DikDjl + DilDjk)
                end
            end
        end
    end
    return W, S, CC
end

#
#
function constitutive(mat::Neohookean, C)
    dim = size(C, 1)
    J2 = det(C)
    if dim == 1
        Jm2n = 1.0 / J2
    elseif dim == 2
        Jm2n = 1.0 / sqrt(J2)
    elseif dim == 3
        Jm2n = 1.0 / cbrt(J2)
    else
        error("unsupported dimension in neohookean: $dim")
    end
    trC = trace(C)
    κ = mat.κ
    μ = mat.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1)
    Wdev = 0.5 * μ * (Jm2n * trC - dim)
    W = Wvol + Wdev
    I = eye(dim)
    IC = C \ I
    Svol = 0.5 * κ * (J2 - 1) * IC
    Sdev = μ * Jm2n * (I - IC * trC / dim)
    S = Svol + Sdev
    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm2n / dim
    CCvol = κ * (J2 * ICxIC - (J2 - 1) * ICoIC)
    CCdev = μJ2n * (trC * (ICxIC / dim + ICoIC) - ox(IC, I) - ox(I, IC))
    CC = CCvol + CCdev
    return W, S, CC
end
