include("minitensor.jl")

using .MiniTensor

struct ScalarField
    name::String
    value::Array{MiniTensor.Scalar}
end

struct VectorField
    name::String
    value::Array{MiniTensor.Vector}
end

struct TensorField
    name::String
    value::Array{MiniTensor.Tensor}
end