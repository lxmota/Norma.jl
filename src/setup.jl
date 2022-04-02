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
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for i ∈ 1 : num_blks
        blk_id = elem_blk_ids[i]
        elem_type = mesh_struct.elem_type(blk_id)
        num_int = default_num_int_pts(elem_type)
        N, dNdξ, w = isoparametric(elem_type, num_int)
    end
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