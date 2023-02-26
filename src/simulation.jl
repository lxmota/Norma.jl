import YAML

include("simulation_def.jl")
include("exodus.jl")
include("model.jl")
include("time.jl")
include("solver.jl")
include("schwarz.jl")


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

function setup_single(params::Dict{Any,Any})
    input_mesh_file = params["input mesh file"]
    output_mesh_file = params["output mesh file"]
    rm(output_mesh_file, force=true)
    output_mesh = Exodus.copy_mesh(input_mesh_file, output_mesh_file)
    params["output_mesh"] = output_mesh
    input_mesh = Exodus.exodus(input_mesh_file)
    params["input_mesh"] = input_mesh
end

function setup_multi(params::Dict{Any,Any})
    domain_names = params["domains"]
    for domain_name ∈ domain_names
        println("domain file: ", domain_name)
        domain_params = setup(domain_name)
        params[domain_name] = domain_params
        domain_params["global_params"] = params
    end
end

function SingleDomainSimulation(params::Dict{Any,Any})
    input_mesh_file = params["input mesh file"]
    output_mesh_file = params["output mesh file"]
    rm(output_mesh_file, force=true)
    output_mesh = Exodus.copy_mesh(input_mesh_file, output_mesh_file)
    params["output_mesh"] = output_mesh
    input_mesh = Exodus.exodus(input_mesh_file)
    params["input_mesh"] = input_mesh
    integrator = create_time_integrator(params)
    solver = create_solver(params)
    model = create_model(params)
    SingleDomainSimulation(params, integrator, solver, model)
end

function MultiDomainSimulation(params::Dict{Any,Any})
    domain_names = params["domains"]
    sub_simulations = Vector{SingleDomainSimulation}()
    for domain_name ∈ domain_names
        println("domain file: ", domain_name)
        domain_params = setup(domain_name)
        params[domain_name] = domain_params
        domain_params["global_params"] = params
        simulation = SingleDomainSimulation(domain_params)
        push!(sub_simulations, simulation)
    end
    MultiDomainSimulation(params, sub_simulations)
end