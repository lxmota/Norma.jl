import YAML

include("exodus.jl")

function setup(input_file::String)
    top_params = YAML.load_file(input_file)
    domain_files = top_params["domains"]
    num_domains = length(domain_files)
    Exodus = exodus_module()
    for domain ∈ 1 : num_domains
        domain_params = YAML.load_file(domain_files[domain])
        domain_mesh = domain_params["mesh"]
        domain_exo = Exodus.exodus(domain_mesh)
        println("domain: ", domain, ", mesh: ", domain_mesh)
    end
end
