module Norma

include("simulation.jl")
include("loop.jl")

function run(input_file::String)
    params = setup(input_file)
    return loop(params)
end

for input_file ∈ ARGS
    run(input_file)
end

end