include("minitensor.jl")

using LinearAlgebra
using .MiniTensor

abstract type Material end
abstract type Solid <: Material end
abstract type Thermal <: Material end

struct SaintVenant_Kirchhoff <: Solid
    E::MTScalar
    ν::MTScalar
    λ::MTScalar
    μ::MTScalar
    ρ::MTScalar
    function SaintVenant_Kirchhoff(E, ν, ρ)
        λ = E * ν / (1.0 + ν) / (1.0 - 2.0 * ν) 
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, λ, μ, ρ)
    end
    function SaintVenant_Kirchhoff(params::Dict{Any, Any})
        E = params["elastic modulus"]
        ν = params["Poisson's ratio"]
        ρ = params["density"]
        new(E, ν, ρ)
    end
end

struct Linear_Elastic <: Solid
    E::MTScalar
    ν::MTScalar
    λ::MTScalar
    μ::MTScalar
    ρ::MTScalar
    function Linear_Elastic(E, ν, ρ)
        λ = E * ν / (1.0 + ν) / (1.0 - 2.0 * ν) 
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, λ, μ, ρ)
    end
    function Linear_Elastic(params::Dict{Any, Any})
        E = params["elastic modulus"]
        ν = params["Poisson's ratio"]
        ρ = params["density"]
        new(E, ν, ρ)
    end
end

struct Neohookean <: Solid
    E::MTScalar
    ν::MTScalar
    κ::MTScalar
    μ::MTScalar
    ρ::MTScalar
    function Neohookean(E, ν, ρ)
        κ = E / (1.0 - 2.0 * ν) / 3.0
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, κ, μ, ρ)
    end
    function Neohookean(params::Dict{Any, Any})
        E = params["elastic modulus"]
        ν = params["Poisson's ratio"]
        ρ = params["density"]
        new(E, ν, ρ)
    end
end

struct Linear_Isotropic <: Thermal
    κ::MTScalar
    function Linear_Isotropic(params::Dict{Any, Any})
        κ = params["thermal conductivity"]
        new(κ)
    end
end

function odot(A, B)
    m = size(A, 1)
    C = zeros(m, m, m, m)
    for a = 1 : m
        for b = 1 : m
            for c = 1 : m
                for d = 1 : m
                    C[a, b, c, d] = A[a, c] * B[b, d] + A[a, d] * B[b, c]
                end
            end
        end
    end
    C = 0.5 * C
    return C
end

function ox(A, B)
    m = size(A, 1)
    C = zeros(m, m, m, m)
    for a = 1 : m
        for b = 1 : m
            for c = 1 : m
                for d = 1 : m
                    C[a, b, c, d] = A[a, b] * B[c, d]
                end
            end
        end
    end
    return C
end

function oxI(A)
    m = size(A, 1)
    C = zeros(m, m, m, m)
    for a = 1 : m
        for b = 1 : m
            for c = 1 : m
                for d = 1 : m
                    C[a, b, c, d] = A[a, b] * I[c, d]
                end
            end
        end
    end
    return C
end

function Iox(B)
    m = size(B, 1)
    C = zeros(m, m, m, m)
    for a = 1 : m
        for b = 1 : m
            for c = 1 : m
                for d = 1 : m
                    C[a, b, c, d] = I[a, b] * B[c, d]
                end
            end
        end
    end
    return C
end

function convect_tangent(CC, S, F)
    n = size(F, 1)
    AA = zeros(n, n, n, n)
    for i = 1 : n
        for j = 1 : n
            for k = 1 : n
                for l = 1 : n
                    s = 0.0
                    for p = 1 : n
                        for q = 1 : n
                            s = s + F[i, p] * CC[p, j, l, q] * F[k, q]
                        end
                    end
                    AA[i, j, k, l] = S[l, j] * I[i, k] + s
                end
            end
        end
    end
    return AA
end

function second_from_fourth(AA)
    dim = size(AA, 1)
    A = zeros(dim * dim)  
    for i = 1 : dim
        for j = 1 : dim
            p = dim * (i - 1) + j
            for k = 1 : dim
                for l = 1 : dim
                    q = dim * (k - 1) + l
                    A[p, q] = AA[i, j, k, l]
                end
            end
        end
    end
    return A
end

function constitutive(material::SaintVenant_Kirchhoff, F::MTTensor)
    dim = size(F, 1)
    C = transpose(F) * F
    E = 0.5 * (C - I)
    λ = material.λ
    μ = material.μ
    trE = tr(E)
    W = 0.5 * λ * trE * trE + μ * tr(E * E)
    S = λ * trE * I + 2.0 * μ * E
    CC = zeros(dim, dim, dim, dim)
    for i = 1 : dim
        for j = 1 : dim
            for k = 1 : dim
                for l = 1 : dim
                    δijδkl = I[i, j] * I[k, l]
                    δikδjl = I[i, k] * I[j, l]
                    δilδjk = I[i, l] * I[j, k]
                    CC[i, j, k, l] = λ * δijδkl + μ * (δikδjl + δilδjk)
                end
            end
        end
    end
    return W, S, CC
end

function constitutive(material::Linear_Elastic, F::MTTensor)
    dim = size(F, 1)
    ∇u = F - I 
    ϵ = symm(∇u)
    λ = material.λ
    μ = material.μ
    trϵ = tr(ϵ)
    W = 0.5 * λ * trϵ * trϵ + μ * tr(ϵ * ϵ)
    σ = λ * trϵ * I + 2.0 * μ * ϵ
    CC = zeros(dim, dim, dim, dim)
    for i = 1 : dim
        for j = 1 : dim
            for k = 1 : dim
                for l = 1 : dim
                    δijδkl = I[i, j] * I[k, l]
                    δikδjl = I[i, k] * I[j, l]
                    δilδjk = I[i, l] * I[j, k]
                    CC[i, j, k, l] = λ * δijδkl + μ * (δikδjl + δilδjk)
                end
            end
        end
    end
    return σ, S, CC
end

function constitutive(material::Neohookean, F::MTTensor)
    dim = size(F, 1)
    C = transpose(F) * F
    J2 = det(C)
    if dim == 1
        Jm2n = 1.0 / J2
    elseif dim == 2
        Jm2n = 1.0 / sqrt(J2)
    elseif dim == 3
        Jm2n = 1.0 / cbrt(J2)
    else
        error("unsupported dimension in neohookean: ", dim)
    end
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1)
    Wdev = 0.5 * μ * (Jm2n * trC - dim)
    W = Wvol + Wdev
    IC = inv(C)
    Svol = 0.5 * κ * (J2 - 1) * IC
    Sdev = μ * Jm2n * (I - IC * trC / dim)
    S = Svol + Sdev
    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm2n / dim
    CCvol = κ * (J2 * ICxIC - (J2 - 1) * ICoIC)
    CCdev = μJ2n * (trC * (ICxIC / dim + ICoIC) - oxI(IC) - Iox(IC))
    CC = CCvol + CCdev
    return W, S, CC
end

function create_material(params::Dict{Any, Any})
    model_name = params["model"]
    if model_name == "linear elastic"
        return Linear_Elastic(params)
    elseif model_name == "Saint-Venant Kirchhoff"
        return SaintVenant_Kirchhoff(params)
    elseif model_name == "neohookean"
        return Neohookean(params)
    else
        error("Unknown material model : ", model_name)
    end
end