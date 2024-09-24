module Norma

include("minitensor.jl")
include("simulation.jl")
include("evolve.jl")

function run(input_file::String)
    sim = create_simulation(input_file)
    evolve(sim)
    return sim
end

function run(params::Dict{Any,Any}, name::String)
    sim = create_simulation(params, name)
    evolve(sim)
    return sim
end

for input_file ∈ ARGS
    run(input_file)
end

end
