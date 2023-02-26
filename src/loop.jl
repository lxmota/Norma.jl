function loop(params::Dict{Any,Any})
    sim_type = params["type"]
    if sim_type == "single"
        return loop_single(params)
    elseif sim_type == "multi"
        return loop_multi(params)
    else
        error("Unknown type of simulation: ", sim_type)
    end
end

function loop_single(params::Dict{Any,Any})
    integrator, solver, model = create_simulation(params)
    apply_ics(model)
    initialize_writing(integrator, model)
    while integrator.time <= integrator.final_time
        model.time = integrator.time
        println("Stop ", integrator.stop, ", time: ", integrator.time)
        apply_bcs(model)
        solve(integrator, model, solver)
        write_step(integrator, model)
        integrator.time = round(integrator.time + integrator.time_step; digits=10)
        integrator.stop += 1
    end
    finalize_writing(model)
    return integrator, solver, model
end

function loop_multi(params::Dict{Any,Any})
    integrators, solvers, models = create_simulations(params)
    apply_ics(models)
    initialize_writing(integrators, models)
    time = initial_time
    stop = 0
    time_step = params["time step"]
    while time <= final_time
        next_time = time + time_step
        synchronize(integrators, models, stop, time)
        println("Stop ", stop, ", time: ", time)
        schwarz(integrators, solvers, models, time, next_time, stop)
        time = round(next_time; digits=10)
        stop += 1
    end
    finalize_writing(models)
    return integrators, solvers, models
end

function schwarz(integrators::Vector{TimeIntegrator}, solvers::Vector{Solver}, models::Vector{Model},
    initial_time::Float64, final_time::Float64, stop::Int64)
    apply_bcs(models)
    if stop == 0
        solve(integrators, solvers, models)
        write_step(integrators, models)
    end
    write_step(integrators, models)
end

function create_simulation(params::Dict{Any,Any})
    integrator = create_time_integrator(params)
    solver = create_solver(params)
    model = create_model(params)
    return integrator, solver, model
end

function create_simulations(params::Dict{Any,Any})
    models = Vector{Model}()
    integrators = Vector{TimeIntegrator}()
    solvers = Vector{Solver}()
    domain_names = params["domains"]
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
        domain_params = params[domain_name]
        integrator_params = domain_params["time integrator"]
        integrator_params["initial time"] = initial_time
        integrator_params["final time"] = final_time
        integrator, solver, model = create_simulation(domain_params)
        model.params["Exodus output interval"] = exodus_interval
        model.params["CSV output interval"] = csv_interval
        push!(models, model)
        push!(solvers, solver)
        push!(integrators, integrator)
    end
    return integrators, solvers, models
end

function apply_ics(models::Vector{Model})
    for model ∈ models
        apply_ics(model)
    end
end

function apply_bcs(models::Vector{Model})
    for model ∈ models
        apply_bcs(model)
    end
end

function solve(integrators::Vector{TimeIntegrator}, solvers::Vector{Solver}, models::Vector{Model})
    num_domains = length(models)
    for domain ∈ 1:num_domains
        solve(integrators[domains], solvers[domain], model[domain])
    end
end

function finalize_writing(models::Vector{Model})
    for model ∈ models
        finalize_writing(model)
    end
end

function initialize_writing(integrators::Vector{TimeIntegrator}, models::Vector{Model})
    num_domains = length(models)
    for domain ∈ 1:num_domains
        initialize_writing(integrators[domains], model[domain])
    end
end

function write_step(integrators::Vector{TimeIntegrator}, models::Vector{Model})
    num_domains = length(models)
    for domain ∈ 1:num_domains
        write_step(integrators[domains], model[domain])
    end
end

function synchronize(integrators::Vector{TimeIntegrator}, models::Vector{Model}, time::Float64, stop::Int64)
    num_domains = length(models)
    for domain ∈ 1:num_domains
        integrators[domain].time = models[domain].time = time
        integrators[domain].stop = stop
    end
end