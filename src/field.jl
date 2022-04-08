include("minitensor.jl")

using .MiniTensor

struct ScalarField
    name::String
    value::Vector{MTScalar}
end

struct VectorField
    name::String
    value::Vector{MTVector}
end

struct TensorField
    name::String
    value::Vector{MTTensor}
end