using LinearAlgebra
using Statistics
using Test

include("../src/norma.jl")
include("helpers.jl")

include("interpolation.jl")
include("single-static-solid-cube.jl")
include("single-static-solid-neumann-bc.jl")
include("single-implicit-dynamic-solid-cube.jl")
include("single-implicit-dynamic-solid-sho.jl")
include("single-implicit-dynamic-solid-clamped.jl")
include("single-explicit-dynamic-solid-cube.jl")
include("single-explicit-dynamic-solid-sho.jl")
include("single-explicit-dynamic-solid-clamped.jl")
include("tet4-static-solid-cube.jl")
include("tet10-static-solid-cube.jl")