include("minitensor.jl")

using .MiniTensor

abstract type Material end

struct SaintVenant_Kirchhoff <: Material
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

struct Linear_Elastic <: Material
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

struct Neohookean <: Material
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

function constitutive(material::SaintVenant_Kirchhoff, F::MTTensor)
    dim = size(F, 1)
    I = eye(dim)
    C = transpose(F) * F
    E = 0.5 * (C - I)
    λ = material.λ
    μ = material.μ
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

function constitutive(material::Linear_Elastic, F::MTTensor)
    dim = size(F, 1)
    I = eye(dim)
    ∇u = F - I 
    E = symm(∇u)
    λ = material.λ
    μ = material.μ
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
    trC = trace(C)
    κ = material.κ
    μ = material.μ
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