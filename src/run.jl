include("setup.jl")
include("loop.jl")

# ARGS = ["cuboids.yaml"]
# cd("/Users/amota/Repos/jlcm/examples/overlap/cuboids")

for input_file ∈ ARGS
    params = setup(input_file)
    loop(params)
end
