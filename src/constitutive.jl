struct SaintVenant_Kirchhoff <: Solid
    E::Float64
    ν::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function SaintVenant_Kirchhoff(params::Dict{Any,Any})
        E = params["elastic modulus"]
        ν = params["Poisson's ratio"]
        ρ = params["density"]
        λ = E * ν / (1.0 + ν) / (1.0 - 2.0 * ν)
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, λ, μ, ρ)
    end
end

struct Linear_Elastic <: Solid
    E::Float64
    ν::Float64
    λ::Float64
    μ::Float64
    ρ::Float64
    function Linear_Elastic(params::Dict{Any,Any})
        E = params["elastic modulus"]
        ν = params["Poisson's ratio"]
        ρ = params["density"]
        λ = E * ν / (1.0 + ν) / (1.0 - 2.0 * ν)
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, λ, μ, ρ)
    end
end

struct Neohookean <: Solid
    E::Float64
    ν::Float64
    κ::Float64
    μ::Float64
    ρ::Float64
    function Neohookean(params::Dict{Any,Any})
        E = params["elastic modulus"]
        ν = params["Poisson's ratio"]
        ρ = params["density"]
        κ = E / (1.0 - 2.0 * ν) / 3.0
        μ = E / (1.0 + ν) / 2.0
        new(E, ν, κ, μ, ρ)
    end
end

struct Linear_Isotropic <: Thermal
    κ::Float64
    function Linear_Isotropic(params::Dict{Any,Any})
        κ = params["thermal conductivity"]
        new(κ)
    end
end

function odot(A::MTTensor, B::MTTensor)
    C = MTTensor4(undef)
    for a ∈ 1:3
        for b ∈ 1:3
            for c ∈ 1:3
                for d ∈ 1:3
                    C[a, b, c, d] = A[a, c] * B[b, d] + A[a, d] * B[b, c]
                end
            end
        end
    end
    C = 0.5 * C
    return C
end

function ox(A::MTTensor, B::MTTensor)
    C = MTTensor4(undef)
    for a ∈ 1:3
        for b ∈ 1:3
            for c ∈ 1:3
                for d ∈ 1:3
                    C[a, b, c, d] = A[a, b] * B[c, d]
                end
            end
        end
    end
    return C
end

function oxI(A::MTTensor)
    C = MTTensor4(undef)
    for a ∈ 1:3
        for b ∈ 1:3
            for c ∈ 1:3
                for d ∈ 1:3
                    C[a, b, c, d] = A[a, b] * I[c, d]
                end
            end
        end
    end
    return C
end

function Iox(B::MTTensor)
    C = MTTensor4(undef)
    for a ∈ 1:3
        for b ∈ 1:3
            for c ∈ 1:3
                for d ∈ 1:3
                    C[a, b, c, d] = I[a, b] * B[c, d]
                end
            end
        end
    end
    return C
end

function convect_tangent(CC::MTTensor4, S::MTTensor, F::MTTensor)
    AA = MTTensor4(undef)
    for i ∈ 1:3
        for j ∈ 1:3
            for k ∈ 1:3
                for l ∈ 1:3
                    s = 0.0
                    for p ∈ 1:3
                        for q ∈ 1:3
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

function second_from_fourth(AA::MTTensor4)
    A = zeros(9, 9)
    for i ∈ 1:3
        for j ∈ 1:3
            p = 3 * (i - 1) + j
            for k ∈ 1:3
                for l ∈ 1:3
                    q = 3 * (k - 1) + l
                    A[p, q] = AA[i, j, k, l]
                end
            end
        end
    end
    return A
end

function constitutive(material::SaintVenant_Kirchhoff, F::MTTensor)
    C = F' * F
    E = 0.5 * (C - I)
    λ = material.λ
    μ = material.μ
    trE = tr(E)
    W = 0.5 * λ * trE * trE + μ * tr(E * E)
    S = λ * trE * I + 2.0 * μ * E
    CC = MTTensor4(undef)
    for i = 1:3
        for j = 1:3
            δᵢⱼ = I[i, j]
            for k = 1:3
                δᵢₖ = I[i, k]; δⱼₖ = I[j, k]
                for l = 1:3
                    δᵢₗ = I[i, l]; δⱼₗ = I[j, l]; δₖₗ = I[k, l]
                    CC[i, j, k, l] = λ * δᵢⱼ * δₖₗ + μ * (δᵢₖ * δⱼₗ + δᵢₗ * δⱼₖ)
                end
            end
        end
    end
    P = F * S
    AA = convect_tangent(CC, S, F)
    return W, P, AA
end

function constitutive(material::Linear_Elastic, F::MTTensor)
    ∇u = F - I
    ϵ = MiniTensor.symm(∇u)
    λ = material.λ
    μ = material.μ
    trϵ = tr(ϵ)
    W = 0.5 * λ * trϵ * trϵ + μ * tr(ϵ * ϵ)
    σ = λ * trϵ * I + 2.0 * μ * ϵ
    CC = MTTensor4(undef)
    for i = 1:3
        for j = 1:3
            δᵢⱼ = I[i, j]
            for k = 1:3
                δᵢₖ = I[i, k]; δⱼₖ = I[j, k]
                for l = 1:3
                    δᵢₗ = I[i, l]; δⱼₗ = I[j, l]; δₖₗ = I[k, l]
                    CC[i, j, k, l] = λ * δᵢⱼ * δₖₗ + μ * (δᵢₖ * δⱼₗ + δᵢₗ * δⱼₖ)
                end
            end
        end
    end
    return W, σ, CC
end

function constitutive(material::Neohookean, F::MTTensor)
    C = F' * F
    J2 = det(C)
    Jm23 = 1.0 / cbrt(J2)
    trC = tr(C)
    κ = material.κ
    μ = material.μ
    Wvol = 0.25 * κ * (J2 - log(J2) - 1)
    Wdev = 0.5 * μ * (Jm23 * trC - 3)
    W = Wvol + Wdev
    IC = inv(C)
    Svol = 0.5 * κ * (J2 - 1) * IC
    Sdev = μ * Jm23 * (I - IC * trC / 3)
    S = Svol + Sdev
    ICxIC = ox(IC, IC)
    ICoIC = odot(IC, IC)
    μJ2n = 2.0 * μ * Jm23 / 3
    CCvol = κ * (J2 * ICxIC - (J2 - 1) * ICoIC)
    CCdev = μJ2n * (trC * (ICxIC / 3 + ICoIC) - oxI(IC) - Iox(IC))
    CC = CCvol + CCdev
    P = F * S
    AA = convect_tangent(CC, S, F)
    return W, P, AA
end

function create_material(params::Dict{Any,Any})
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

function get_p_wave_modulus(material::Solid)
    return material.λ + 2.0 * material.μ
end