module Norma

include("simulation.jl")
include("loop.jl")

function run(input_file::String)
    simulation = create_simulation(input_file)
    loop(simulation)
    return simulation
end

for input_file ∈ ARGS
    run(input_file)
end

end