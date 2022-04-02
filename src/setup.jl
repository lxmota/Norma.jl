import YAML

include("exodus.jl")
include("interpolation.jl")

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
    elem_blk_ids = mesh_struct.get_elem_blk_ids()
    num_blks = length(elem_blk_ids)
    for i ∈ 1 : num_blks
        blk_id = elem_blk_ids[i]
        elem_type = mesh_struct.elem_type(blk_id)
        num_int = default_num_int_pts(elem_type)
        N, dNdξ, w = isoparametric(elem_type, num_int)
    end
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
