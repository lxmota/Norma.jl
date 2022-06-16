import YAML

include("exodus.jl")
include("interpolation.jl")

Exodus = exodus_module()

function setup(input_file::String)
    params = YAML.load_file(input_file)
    sim_type = params["type"]
    if sim_type == "single"
        setup_single(params)
    elseif sim_type == "multi"
        setup_multi(params)
    else
        error("Unknown type of simulation: ", sim_type)
    end
    return params
end

function setup_single(params::Dict{Any, Any})
    mesh_file = params["mesh"]
    mesh_struct = Exodus.exodus(mesh_file)
    params["mesh_struct"] = mesh_struct
end

function setup_multi(params::Dict{Any, Any})
    domain_files = params["domains"]
    num_domains = length(domain_files)
    for domain ∈ 1 : num_domains
        domain_file = domain_files[domain]
        println("domain: ", domain, ", domain file: ", domain_file)
        setup(domain_file)
    end
end