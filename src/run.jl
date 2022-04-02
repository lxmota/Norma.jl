include("setup.jl")

# ARGS = ["cuboids.yaml"]
# cd("/Users/amota/Repos/jlcm/examples/overlap/cuboids")

for input_file ∈ ARGS
    setup(input_file)
end
