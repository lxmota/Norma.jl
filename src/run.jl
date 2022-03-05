import YAML

function setup(input_file::String)
    top_params = YAML.load_file(input_file)
    domain_files = top_params["domains"]
    num_domains = length(domain_files)
    println(num_domains)
end
