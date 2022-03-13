import YAML

include("exodus.jl")

Exodus = exodus_module()

function setup(input_file::String)
    params = YAML.load_file(input_file)
    sim_type = params["type"]
    if sim_type == "single"
        setup_single(input_file)
    elseif sim_type == "multi"
        setup_multi(input_file)
    else
        error("Unknown type of simulation: ", sim_type)
    end
end

function setup_single(input_file::String)
    params = YAML.load_file(input_file)
    mesh_file = params["mesh"]
    mesh_struct = Exodus.exodus(mesh_file)
end

function setup_multi(input_file::String)
    params = YAML.load_file(input_file)
    domain_files = params["domains"]
    num_domains = length(domain_files)
    for domain ∈ 1 : num_domains
        println("domain: ", domain)
        domain_file = domain_files[domain]
        setup(domain_file)
    end
end

for input_file ∈ ARGS
    setup(input_file)
end
