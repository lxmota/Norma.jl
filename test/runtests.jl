using Statistics
using Test

include("../src/norma.jl")
include("helpers.jl")

@testset "Norma" begin
    include("single-static-solid-cube.jl")
    include("single-dynamic-solid-cube.jl")
end