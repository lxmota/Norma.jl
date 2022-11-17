module JLCM

include("setup.jl")
include("init.jl")
include("loop.jl")

for input_file ∈ ARGS
    params = setup(input_file)
    init(params)
    loop(params)
end

end