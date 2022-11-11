abstract type Solution end
abstract type QuasiStatic <: Solution end
abstract type TransientDynamic <: Solution end

mutable struct Continuation <: QuasiStatic
end