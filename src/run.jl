module JLCM

include("setup.jl")
include("init.jl")
include("loop.jl")

ARGS = ["cube.yaml"]
# cd("/Users/amota/Repos/jlcm/examples/single/static-solid/unit-cube")

for input_file ∈ ARGS
    params = setup(input_file)
    init(params)
    loop(params)
end

end