module Norma

include("setup.jl")
include("init.jl")
include("loop.jl")

function run(input_file::String)
    params = setup(input_file)
    init(params)
    integrator, solver, model = loop(params)
    return integrator, solver, model
end

for input_file ∈ ARGS
    run(input_file)
end

end