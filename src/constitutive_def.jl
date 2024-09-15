using Enzyme
using LinearAlgebra

include("minitensor.jl")

using .MiniTensor

abstract type Material end
abstract type Solid <: Material end
abstract type Thermal <: Material end
