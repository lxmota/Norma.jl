module Norma

include("simulation.jl")
include("evolve.jl")

function run(input_file::String)
    simulation = create_simulation(input_file)
    evolve(simulation)
    return simulation
end

for input_file ∈ ARGS
    run(input_file)
end

end