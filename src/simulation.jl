import YAML

include("simulation_def.jl")
include("exodus.jl")
include("model.jl")
include("time.jl")
include("solver.jl")
include("schwarz.jl")


Exodus = exodus_module()
Exodus.SHOW_BANNER = false

function create_simulation(input_file::String)
    println("Reading simulation file: ", input_file)
    params = YAML.load_file(input_file)
    params["name"] = input_file
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
    name = params["name"]
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
    SingleDomainSimulation(name, params, integrator, solver, model)
end

function MultiDomainSimulation(params::Dict{Any,Any})
    name = params["name"]
    domain_names = params["domains"]
    subsims = Vector{SingleDomainSimulation}()
    initial_time = params["initial time"]
    final_time = params["final time"]
    exodus_interval = get(params, "Exodus output interval", 1)
    csv_interval = get(params, "CSV output interval", 0)
    sim_type = "none"
    for domain_name ∈ domain_names
        println("Reading subsimulation file: ", domain_name)
        subparams = YAML.load_file(domain_name)
        subparams["name"] = domain_name
        subparams["time integrator"]["initial time"] = initial_time
        subparams["time integrator"]["final time"] = final_time
        subparams["Exodus output interval"] = exodus_interval
        subparams["CSV output interval"] = csv_interval
        subsim = SingleDomainSimulation(subparams)
        params[domain_name] = subsim.params
        subsim.params["global_params"] = params
        integrator_name = subsim.params["time integrator"]["type"]
        subsim_type = is_static_or_dynamic(integrator_name)
        if sim_type == "none"
            sim_type = subsim_type
        elseif subsim_type ≠ sim_type
            error("Multidomain subdomains must be all static or dynamic")
        end
        params["subdomains type"] = sim_type
        push!(subsims, subsim)
    end
    schwarz_controller = create_schwarz_controller(params)
    MultiDomainSimulation(name, params, schwarz_controller, subsims)
end