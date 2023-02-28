import YAML

include("simulation_def.jl")
include("exodus.jl")
include("model.jl")
include("time.jl")
include("solver.jl")
include("schwarz.jl")


Exodus = exodus_module()

function create_simulation(input_file::String)
    params = YAML.load_file(input_file)
    sim_type = params["type"]
    if sim_type == "single"
        return SingleDomainSimulation(params)
    elseif sim_type == "multi"
        return MultiDomainSimulation(params)
    else
        error("Unknown type of simulation: ", sim_type)
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
    initial_time = params["initial time"]
    final_time = params["final time"]
    exodus_interval = 1
    if haskey(params, "Exodus output interval") == true
        exodus_interval = params["Exodus output interval"]
    end
    csv_interval = 0
    if haskey(params, "CSV output interval") == true
        csv_interval = params["CSV output interval"]
    end
    for domain_name ∈ domain_names
        println("domain file: ", domain_name)
        domain_params = setup(domain_name)
        params[domain_name] = domain_params
        domain_params["global_params"] = params
        integrator_params = domain_params["time integrator"]
        integrator_params["initial time"] = initial_time
        integrator_params["final time"] = final_time
        simulation = SingleDomainSimulation(domain_params)
        simulation.params["Exodus output interval"] = exodus_interval
        simulation.params["CSV output interval"] = csv_interval
        push!(sub_simulations, simulation)
    end
    MultiDomainSimulation(params, sub_simulations)
end