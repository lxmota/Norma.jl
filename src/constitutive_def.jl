using LinearAlgebra

using .MiniTensor

abstract type Material end
abstract type Solid <: Material end
abstract type Thermal <: Material end
abstract type OperatorInference <: Material end
abstract type ConvectionDiffusion <: Material end
